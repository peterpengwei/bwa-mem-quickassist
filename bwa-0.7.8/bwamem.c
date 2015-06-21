#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

// [QA] For HLS
#include <fcntl.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
// [QA] For HLS

#include "kstring.h"
#include "bwamem.h"
#include "bntseq.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

// [QA] Macro for memory space size (Unit: Byte)
#define INPUT_SPACE_SIZE (1<<18)
#define INPUT_THRESHOLD 1024
#define OUTPUT_SPACE_SIZE (1<<14)

#define MIN_BATCH_SIZE 64

#include "batch.h"

batch* reqTBBSpace(); 
void releaseBatchSpace(batch*);

/* Theory on probability and scoring *ungapped* alignment
 *
 * s'(a,b) = log[P(b|a)/P(b)] = log[4P(b|a)], assuming uniform base distribution
 * s'(a,a) = log(4), s'(a,b) = log(4e/3), where e is the error rate
 *
 * Scale s'(a,b) to s(a,a) s.t. s(a,a)=x. Then s(a,b) = x*s'(a,b)/log(4), or conversely: s'(a,b)=s(a,b)*log(4)/x
 *
 * If the matching score is x and mismatch penalty is -y, we can compute error rate e:
 *   e = .75 * exp[-log(4) * y/x]
 *
 * log P(seq) = \sum_i log P(b_i|a_i) = \sum_i {s'(a,b) - log(4)}
 *   = \sum_i { s(a,b)*log(4)/x - log(4) } = log(4) * (S/x - l)
 *
 * where S=\sum_i s(a,b) is the alignment score. Converting to the phred scale:
 *   Q(seq) = -10/log(10) * log P(seq) = 10*log(4)/log(10) * (l - S/x) = 6.02 * (l - S/x)
 *
 *
 * Gap open (zero gap): q' = log[P(gap-open)], r' = log[P(gap-ext)] (see Durbin et al. (1998) Section 4.1)
 * Then q = x*log[P(gap-open)]/log(4), r = x*log[P(gap-ext)]/log(4)
 *
 * When there are gaps, l should be the length of alignment matches (i.e. the M operator in CIGAR)
 */

mem_opt_t *mem_opt_init()
{
	mem_opt_t *o;
	o = calloc(1, sizeof(mem_opt_t));
	o->flag = 0;
	o->a = 1; o->b = 4;
	o->o_del = o->o_ins = 6;
	o->e_del = o->e_ins = 1;
	o->w = 100;
	o->T = 30;
	o->zdrop = 100;
	o->pen_unpaired = 17;
	o->pen_clip5 = o->pen_clip3 = 5;
	o->min_seed_len = 19;
	o->split_width = 10;
	o->max_occ = 10000;
	o->max_chain_gap = 10000;
	o->max_ins = 10000;
	o->mask_level = 0.50;
	o->chain_drop_ratio = 0.50;
	o->split_factor = 1.5;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->batch_size = 1; // [QA] initialize batch size
	o->max_matesw = 100;
	o->mask_level_redun = 0.95;
	o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
//	o->mapQ_coef_len = o->mapQ_coef_fac = 0;
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}

/***************************
 * SMEM iterator interface *
 ***************************/

struct __smem_i {
	const bwt_t *bwt;
	const uint8_t *query;
	int start, len;
	bwtintv_v *matches; // matches; to be returned by smem_next()
	bwtintv_v *sub;     // sub-matches inside the longest match; temporary
	bwtintv_v *tmpvec[2]; // temporary arrays
};

smem_i *smem_itr_init(const bwt_t *bwt)
{
	smem_i *itr;
	itr = calloc(1, sizeof(smem_i));
	itr->bwt = bwt;
	itr->tmpvec[0] = calloc(1, sizeof(bwtintv_v));
	itr->tmpvec[1] = calloc(1, sizeof(bwtintv_v));
	itr->matches   = calloc(1, sizeof(bwtintv_v));
	itr->sub       = calloc(1, sizeof(bwtintv_v));
	return itr;
}

void smem_itr_destroy(smem_i *itr)
{
	free(itr->tmpvec[0]->a); free(itr->tmpvec[0]);
	free(itr->tmpvec[1]->a); free(itr->tmpvec[1]);
	free(itr->matches->a);   free(itr->matches);
	free(itr->sub->a);       free(itr->sub);
	free(itr);
}

void smem_set_query(smem_i *itr, int len, const uint8_t *query)
{
	itr->query = query;
	itr->start = 0;
	itr->len = len;
}

const bwtintv_v *smem_next2(smem_i *itr, int split_len, int split_width, int start_width)
{
	int i, max, max_i, ori_start;
	itr->tmpvec[0]->n = itr->tmpvec[1]->n = itr->matches->n = itr->sub->n = 0;
	if (itr->start >= itr->len || itr->start < 0) return 0;
	while (itr->start < itr->len && itr->query[itr->start] > 3) ++itr->start; // skip ambiguous bases
	if (itr->start == itr->len) return 0;
	ori_start = itr->start;
	itr->start = bwt_smem1(itr->bwt, itr->len, itr->query, ori_start, start_width, itr->matches, itr->tmpvec); // search for SMEM
	if (itr->matches->n == 0) return itr->matches; // well, in theory, we should never come here
	for (i = max = 0, max_i = 0; i < itr->matches->n; ++i) { // look for the longest match
		bwtintv_t *p = &itr->matches->a[i];
		int len = (uint32_t)p->info - (p->info>>32);
		if (max < len) max = len, max_i = i;
	}
	if (split_len > 0 && max >= split_len && itr->matches->a[max_i].x[2] <= split_width) { // if the longest SMEM is unique and long
		int j;
		bwtintv_v *a = itr->tmpvec[0]; // reuse tmpvec[0] for merging
		bwtintv_t *p = &itr->matches->a[max_i];
		bwt_smem1(itr->bwt, itr->len, itr->query, ((uint32_t)p->info + (p->info>>32))>>1, itr->matches->a[max_i].x[2]+1, itr->sub, itr->tmpvec); // starting from the middle of the longest MEM
		i = j = 0; a->n = 0;
		while (i < itr->matches->n && j < itr->sub->n) { // ordered merge
			int64_t xi = itr->matches->a[i].info>>32<<32 | (itr->len - (uint32_t)itr->matches->a[i].info);
			int64_t xj = itr->sub->a[j].info>>32<<32 | (itr->len - (uint32_t)itr->sub->a[j].info);
			if (xi < xj) {
				kv_push(bwtintv_t, *a, itr->matches->a[i]);
				++i;
			} else if ((uint32_t)itr->sub->a[j].info - (itr->sub->a[j].info>>32) >= max>>1 && (uint32_t)itr->sub->a[j].info > ori_start) {
				kv_push(bwtintv_t, *a, itr->sub->a[j]);
				++j;
			} else ++j;
		}
		for (; i < itr->matches->n; ++i) kv_push(bwtintv_t, *a, itr->matches->a[i]);
		for (; j < itr->sub->n; ++j)
			if ((uint32_t)itr->sub->a[j].info - (itr->sub->a[j].info>>32) >= max>>1 && (uint32_t)itr->sub->a[j].info > ori_start)
				kv_push(bwtintv_t, *a, itr->sub->a[j]);
		kv_copy(bwtintv_t, *itr->matches, *a);
	}
	return itr->matches;
}

const bwtintv_v *smem_next(smem_i *itr, int split_len, int split_width)
{
	return smem_next2(itr, split_len, split_width, 1);
}

/********************************
 * Chaining while finding SMEMs *
 ********************************/

typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
} mem_seed_t;

typedef struct {
	int n, m;
	int64_t pos;
	mem_seed_t *seeds;
} mem_chain_t;

typedef struct { size_t n, m; mem_chain_t *a;  } mem_chain_v;

#include "kbtree.h"

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, mem_chain_t, chain_cmp)

static int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c, const mem_seed_t *p)
{
	int64_t qend, rend, x, y;
	const mem_seed_t *last = &c->seeds[c->n-1];
	qend = last->qbeg + last->len;
	rend = last->rbeg + last->len;
	if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend && p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
		return 1; // contained seed; do nothing
	if ((last->rbeg < l_pac || c->seeds[0].rbeg < l_pac) && p->rbeg >= l_pac) return 0; // don't chain if on different strand
	x = p->qbeg - last->qbeg; // always non-negtive
	y = p->rbeg - last->rbeg;
	if (y >= 0 && x - y <= opt->w && y - x <= opt->w && x - last->len < opt->max_chain_gap && y - last->len < opt->max_chain_gap) { // grow the chain
		if (c->n == c->m) {
			c->m <<= 1;
			c->seeds = realloc(c->seeds, c->m * sizeof(mem_seed_t));
		}
		c->seeds[c->n++] = *p;
		return 1;
	}
	return 0; // request to add a new chain
}

static void mem_insert_seed(const mem_opt_t *opt, int64_t l_pac, kbtree_t(chn) *tree, smem_i *itr)
{
	const bwtintv_v *a;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	int start_width = (opt->flag & MEM_F_NO_EXACT)? 2 : 1;
	split_len = split_len < itr->len? split_len : itr->len;
	while ((a = smem_next2(itr, split_len, opt->split_width, start_width)) != 0) { // to find all SMEM and some internal MEM
		int i;
		for (i = 0; i < a->n; ++i) { // go through each SMEM/MEM up to itr->start
			bwtintv_t *p = &a->a[i];
			int slen = (uint32_t)p->info - (p->info>>32); // seed length
			int64_t k;
			if (slen < opt->min_seed_len || p->x[2] > opt->max_occ) continue; // ignore if too short or too repetitive
			for (k = 0; k < p->x[2]; ++k) {
				mem_chain_t tmp, *lower, *upper;
				mem_seed_t s;
				int to_add = 0;
				s.rbeg = tmp.pos = bwt_sa(itr->bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
				s.qbeg = p->info>>32;
				s.len  = slen;
				if (bwa_verbose >= 5) printf("* Found SEED: length=%d,query_beg=%d,ref_beg=%ld\n", s.len, s.qbeg, (long)s.rbeg);
				if (s.rbeg < l_pac && l_pac < s.rbeg + s.len) continue; // bridging forward-reverse boundary; skip
				if (kb_size(tree)) {
					kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
					if (!lower || !test_and_merge(opt, l_pac, lower, &s)) to_add = 1;
				} else to_add = 1;
				if (to_add) { // add the seed as a new chain
					tmp.n = 1; tmp.m = 4;
					tmp.seeds = calloc(tmp.m, sizeof(mem_seed_t));
					tmp.seeds[0] = s;
					kb_putp(chn, tree, &tmp);
				}
			}
		}
	}
}

int mem_chain_weight(const mem_chain_t *c)
{
	int64_t end;
	int j, w = 0, tmp;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->qbeg >= end) w += s->len;
		else if (s->qbeg + s->len > end) w += s->qbeg + s->len - end;
		end = end > s->qbeg + s->len? end : s->qbeg + s->len;
	}
	tmp = w;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->rbeg >= end) w += s->len;
		else if (s->rbeg + s->len > end) w += s->rbeg + s->len - end;
		end = end > s->qbeg + s->len? end : s->qbeg + s->len;
	}
	return w < tmp? w : tmp;
}

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn)
{
	int i, j;
	for (i = 0; i < chn->n; ++i) {
		mem_chain_t *p = &chn->a[i];
		err_printf("* Found CHAIN(%d): n=%d; weight=%d", i, p->n, mem_chain_weight(p));
		for (j = 0; j < p->n; ++j) {
			bwtint_t pos;
			int is_rev, ref_id;
			pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
			if (is_rev) pos -= p->seeds[j].len - 1;
			bns_cnt_ambi(bns, pos, p->seeds[j].len, &ref_id);
			err_printf("\t%d;%d,%ld(%s:%c%ld)", p->seeds[j].len, p->seeds[j].qbeg, (long)p->seeds[j].rbeg, bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - bns->anns[ref_id].offset) + 1);
		}
		err_putchar('\n');
	}
}

mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, int64_t l_pac, int len, const uint8_t *seq)
{
	mem_chain_v chain;
	smem_i *itr;
	kbtree_t(chn) *tree;

	kv_init(chain);
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	tree = kb_init(chn, KB_DEFAULT_SIZE);
	itr = smem_itr_init(bwt);
	smem_set_query(itr, len, seq);
	mem_insert_seed(opt, l_pac, tree, itr);

	kv_resize(mem_chain_t, chain, kb_size(tree));

	#define traverse_func(p_) (chain.a[chain.n++] = *(p_))
	__kb_traverse(mem_chain_t, tree, traverse_func);
	#undef traverse_func

	smem_itr_destroy(itr);
	kb_destroy(chn, tree);
	return chain;
}

/********************
 * Filtering chains *
 ********************/

typedef struct {
	int beg, end, w;
	void *p, *p2;
} flt_aux_t;

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(mem_flt, flt_aux_t, flt_lt)

int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *chains)
{
	flt_aux_t *a;
	int i, j, n;
	if (n_chn <= 1) return n_chn; // no need to filter
	a = malloc(sizeof(flt_aux_t) * n_chn);
	for (i = 0; i < n_chn; ++i) {
		mem_chain_t *c = &chains[i];
		int w;
		w = mem_chain_weight(c);
		a[i].beg = c->seeds[0].qbeg;
		a[i].end = c->seeds[c->n-1].qbeg + c->seeds[c->n-1].len;
		a[i].w = w; a[i].p = c; a[i].p2 = 0;
	}
	ks_introsort(mem_flt, n_chn, a);
	{ // reorder chains such that the best chain appears first
		mem_chain_t *swap;
		swap = malloc(sizeof(mem_chain_t) * n_chn);
		for (i = 0; i < n_chn; ++i) {
			swap[i] = *((mem_chain_t*)a[i].p);
			a[i].p = &chains[i]; // as we will memcpy() below, a[i].p is changed
		}
		memcpy(chains, swap, sizeof(mem_chain_t) * n_chn);
		free(swap);
	}
	for (i = 1, n = 1; i < n_chn; ++i) {
		for (j = 0; j < n; ++j) {
			int b_max = a[j].beg > a[i].beg? a[j].beg : a[i].beg;
			int e_min = a[j].end < a[i].end? a[j].end : a[i].end;
			if (e_min > b_max) { // have overlap
				int min_l = a[i].end - a[i].beg < a[j].end - a[j].beg? a[i].end - a[i].beg : a[j].end - a[j].beg;
				if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
					if (a[j].p2 == 0) a[j].p2 = a[i].p;
					if (a[i].w < a[j].w * opt->chain_drop_ratio && a[j].w - a[i].w >= opt->min_seed_len<<1)
						break;
				}
			}
		}
		if (j == n) a[n++] = a[i]; // if have no significant overlap with better chains, keep it.
	}
	for (i = 0; i < n; ++i) { // mark chains to be kept
		mem_chain_t *c = (mem_chain_t*)a[i].p;
		if (c->n > 0) c->n = -c->n;
		c = (mem_chain_t*)a[i].p2;
		if (c && c->n > 0) c->n = -c->n;
	}
	free(a);
	for (i = 0; i < n_chn; ++i) { // free discarded chains
		mem_chain_t *c = &chains[i];
		if (c->n >= 0) {
			free(c->seeds);
			c->n = c->m = 0;
		} else c->n = -c->n;
	}
	for (i = n = 0; i < n_chn; ++i) { // squeeze out discarded chains
		if (chains[i].n > 0) {
			if (n != i) chains[n++] = chains[i];
			else ++n;
		}
	}
	return n;
}

/******************************
 * De-overlap single-end hits *
 ******************************/

#define alnreg_slt2(a, b) ((a).re < (b).re)
KSORT_INIT(mem_ars2, mem_alnreg_t, alnreg_slt2)

#define alnreg_slt(a, b) ((a).score > (b).score || ((a).score == (b).score && ((a).rb < (b).rb || ((a).rb == (b).rb && (a).qb < (b).qb))))
KSORT_INIT(mem_ars, mem_alnreg_t, alnreg_slt)

#define alnreg_hlt(a, b) ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))
KSORT_INIT(mem_ars_hash, mem_alnreg_t, alnreg_hlt)

int mem_sort_and_dedup(int n, mem_alnreg_t *a, float mask_level_redun)
{
	int m, i, j;
	if (n <= 1) return n;
	ks_introsort(mem_ars2, n, a);
	for (i = 1; i < n; ++i) {
		mem_alnreg_t *p = &a[i];
		if (p->rb >= a[i-1].re) continue;
		for (j = i - 1; j >= 0 && p->rb < a[j].re; --j) {
			mem_alnreg_t *q = &a[j];
			int64_t or, oq, mr, mq;
			if (q->qe == q->qb) continue; // a[j] has been excluded
			or = q->re - p->rb; // overlap length on the reference
			oq = q->qb < p->qb? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
			mr = q->re - q->rb < p->re - p->rb? q->re - q->rb : p->re - p->rb; // min ref len in alignment
			mq = q->qe - q->qb < p->qe - p->qb? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
			if (or > mask_level_redun * mr && oq > mask_level_redun * mq) { // one of the hits is redundant
				if (p->score < q->score) {
					p->qe = p->qb;
					break;
				} else q->qe = q->qb;
			}
		}
	}
	for (i = 0, m = 0; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	n = m;
	ks_introsort(mem_ars, n, a);
	for (i = 1; i < n; ++i) { // mark identical hits
		if (a[i].score == a[i-1].score && a[i].rb == a[i-1].rb && a[i].qb == a[i-1].qb)
			a[i].qe = a[i].qb;
	}
	for (i = 1, m = 1; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	return m;
}

int mem_test_and_remove_exact(const mem_opt_t *opt, int n, mem_alnreg_t *a, int qlen)
{
	if (!(opt->flag & MEM_F_NO_EXACT) || n == 0 || a->truesc != qlen * opt->a) return n;
	memmove(a, a + 1, (n - 1) * sizeof(mem_alnreg_t));
	return n - 1;
}

void mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id) // IMPORTANT: must run mem_sort_and_dedup() before calling this function
{ // similar to the loop in mem_chain_flt()
	int i, k, tmp;
	kvec_t(int) z;
	if (n == 0) return;
	kv_init(z);
	for (i = 0; i < n; ++i) a[i].sub = 0, a[i].secondary = -1, a[i].hash = hash_64(id+i);
	ks_introsort(mem_ars_hash, n, a);
	tmp = opt->a + opt->b;
	tmp = opt->o_del + opt->e_del > tmp? opt->o_del + opt->e_del : tmp;
	tmp = opt->o_ins + opt->e_ins > tmp? opt->o_ins + opt->e_ins : tmp;
	kv_push(int, z, 0);
	for (i = 1; i < n; ++i) {
		for (k = 0; k < z.n; ++k) {
			int j = z.a[k];
			int b_max = a[j].qb > a[i].qb? a[j].qb : a[i].qb;
			int e_min = a[j].qe < a[i].qe? a[j].qe : a[i].qe;
			if (e_min > b_max) { // have overlap
				int min_l = a[i].qe - a[i].qb < a[j].qe - a[j].qb? a[i].qe - a[i].qb : a[j].qe - a[j].qb;
				if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
					if (a[j].sub == 0) a[j].sub = a[i].score;
					if (a[j].score - a[i].score <= tmp) ++a[j].sub_n;
					break;
				}
			}
		}
		if (k == z.n) kv_push(int, z, i);
		else a[i].secondary = z.a[k];
	}
	free(z.a);
}

/****************************************
 * Construct the alignment from a chain *
 ****************************************/

/* mem_chain2aln() vs mem_chain2aln_short()
 *
 * mem_chain2aln() covers all the functionality of mem_chain2aln_short().
 * However, it may waste time on extracting the reference sequences given a
 * very long query. mem_chain2aln_short() is faster for very short chains in a
 * long query. It may fail when the matches are long or reach the end of the
 * query. In this case, mem_chain2aln() will be called again.
 * mem_chain2aln_short() is almost never used for short-read alignment.
 */

#define MEM_SHORT_EXT 50
#define MEM_SHORT_LEN 200
#define MAX_BAND_TRY  2

int mem_chain2aln_short(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av)
{
	int i, qb, qe, xtra;
	int64_t rb, re, rlen;
	uint8_t *rseq = 0;
	mem_alnreg_t a;
	kswr_t x;

	if (c->n == 0) return -1;
	qb = l_query;  qe = 0;
	rb = l_pac<<1; re = 0;
	memset(&a, 0, sizeof(mem_alnreg_t));
	for (i = 0; i < c->n; ++i) {
		const mem_seed_t *s = &c->seeds[i];
		qb = qb < s->qbeg? qb : s->qbeg;
		qe = qe > s->qbeg + s->len? qe : s->qbeg + s->len;
		rb = rb < s->rbeg? rb : s->rbeg;
		re = re > s->rbeg + s->len? re : s->rbeg + s->len;
		a.seedcov += s->len;
	}
	qb -= MEM_SHORT_EXT; qe += MEM_SHORT_EXT;
	if (qb <= 10 || qe >= l_query - 10) return 1; // because ksw_align() does not support end-to-end alignment
	rb -= MEM_SHORT_EXT; re += MEM_SHORT_EXT;
	rb = rb > 0? rb : 0;
	re = re < l_pac<<1? re : l_pac<<1;
	if (rb < l_pac && l_pac < re) {
		if (c->seeds[0].rbeg < l_pac) re = l_pac;
		else rb = l_pac;
	}
	if ((re - rb) - (qe - qb) > MEM_SHORT_EXT || (qe - qb) - (re - rb) > MEM_SHORT_EXT) return 1;
	if (qe - qb >= opt->w * 4 || re - rb >= opt->w * 4) return 1;
	if (qe - qb >= MEM_SHORT_LEN || re - rb >= MEM_SHORT_LEN) return 1;

	rseq = bns_get_seq(l_pac, pac, rb, re, &rlen);
	assert(rlen == re - rb);
	xtra = KSW_XSUBO | KSW_XSTART | ((qe - qb) * opt->a < 250? KSW_XBYTE : 0) | (opt->min_seed_len * opt->a);
	x = ksw_align2(qe - qb, (uint8_t*)query + qb, re - rb, rseq, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, xtra, 0);
	free(rseq);
	if (x.tb < MEM_SHORT_EXT>>1 || x.te > re - rb - (MEM_SHORT_EXT>>1)) return 1;

	a.rb = rb + x.tb; a.re = rb + x.te + 1;
	a.qb = qb + x.qb; a.qe = qb + x.qe + 1;
	a.score = x.score;
	a.csub = x.score2;
	kv_push(mem_alnreg_t, *av, a);
	if (bwa_verbose >= 4) printf("** Added alignment region via mem_chain2aln_short(): [%d,%d) <=> [%ld,%ld)\n", a.qb, a.qe, (long)a.rb, (long)a.re);
	return 0;
}

static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

typedef struct {
	uint8_t* leftQs;
	int16_t leftQlen;
	uint8_t* rightQs;
	int16_t rightQlen;
	uint8_t* leftRs;
	int16_t leftRlen;
	uint8_t* rightRs;
	int16_t rightRlen;
	int16_t h0;
	int16_t regScore;
	int16_t qBeg; 	
	int32_t idx;
	mem_alnreg_t* reg;
	int valid;
	int string_size;
	int32_t qe_offset;
	int64_t rbeg_offset;
	int64_t re_offset;
} ext_param_t;

typedef struct {
	int32_t idx;
	int16_t qBeg;
	int16_t rBeg;
	int16_t qEnd;
	int16_t rEnd;
	int16_t score;
	int16_t trueScore;
	int16_t width;
} ext_res_t;

typedef struct {
	int64_t rmax[2];
	uint64_t* srt;
	uint8_t* rseq;
} sw_pre_result_t;

// [QA] the batched-processing version of mem_chain2aln
void mem_chain2aln_preprocess(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, sw_pre_result_t* pre_result)
{
	int i;
	//int i, max_off[2];
	int64_t rlen, rmax[2], max = 0;
	uint8_t *rseq = 0;
	uint64_t *srt;

	assert(c->n > 0);
	if (c->n == 0) return;
	// get the max possible span
	rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) {
		int64_t b, e;
		const mem_seed_t *t = &c->seeds[i];
		b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
		e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
		rmax[0] = rmax[0] < b? rmax[0] : b;
		rmax[1] = rmax[1] > e? rmax[1] : e;
		if (t->len > max) max = t->len;
	}
	rmax[0] = rmax[0] > 0? rmax[0] : 0;
	rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
	if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
		if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
		else rmax[0] = l_pac;
	}
	// retrieve the reference sequence
	rseq = bns_get_seq(l_pac, pac, rmax[0], rmax[1], &rlen);
	assert(rlen == rmax[1] - rmax[0]);

	srt = malloc(c->n * 8);
	for (i = 0; i < c->n; ++i)
		srt[i] = (uint64_t)c->seeds[i].len<<32 | i;
	ks_introsort_64(c->n, srt);
	free (pre_result->srt);
	free (pre_result->rseq);
	pre_result->srt = srt;
	pre_result->rseq = rseq;
	pre_result->rmax[0] = rmax[0];
	pre_result->rmax[1] = rmax[1];
}

void original_extension(mem_alnreg_t* reg_dut, ext_param_t* cur_param, const mem_opt_t* opt) {
	mem_alnreg_t* a = reg_dut;
	int i;
	int aw[2];
	aw[0] = aw[1] = opt->w;
	int max_off[2];
	if (cur_param->leftQlen > 0) { // left extension
		int qle, tle, gtle, gscore;
		for (i = 0; i < MAX_BAND_TRY; ++i) {
			int prev = a->score;
			aw[0] = opt->w << i;
			a->score = ksw_extend2(cur_param->leftQlen, cur_param->leftQs, cur_param->leftRlen, cur_param->leftRs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, cur_param->h0, &qle, &tle, &gtle, &gscore, &max_off[0]);
			if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
		}
		if (bwa_verbose >= 4) printf("left_ext: qle=%d, tle=%d, gtle=%d, gscore=%d, max_off[0]=%d, score=%d\n", qle, tle, gtle, gscore, max_off[0], a->score);
		if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
			a->qb = cur_param->qBeg - qle, a->rb = cur_param->rbeg_offset - tle;
			a->truesc = a->score;
		} else { // to-end extension
			a->qb = 0, a->rb = cur_param->rbeg_offset - gtle;
			a->truesc = gscore;
		}
	} 

	if (cur_param->rightQlen > 0) { // right extension
		int qle, tle, gtle, gscore, sc0 = a->score;
		for (i = 0; i < MAX_BAND_TRY; ++i) {
			int prev = a->score;
			aw[1] = opt->w << i;
			a->score = ksw_extend2(cur_param->rightQlen, cur_param->rightQs, cur_param->rightRlen, cur_param->rightRs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
			if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
		}
		if (bwa_verbose >= 4) printf("right_ext: qle=%d, tle=%d, gtle=%d, gscore=%d, max_off[1]=%d, score=%d\n", qle, tle, gtle, gscore, max_off[1], a->score);
		// similar to the above
		if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
			a->qe = cur_param->qe_offset + qle, a->re = cur_param->re_offset + tle;
			a->truesc += a->score - sc0;
		} else { // to-end extension
			a->qe = cur_param->qe_offset + cur_param->rightQlen, a->re = cur_param->re_offset + gtle;
			a->truesc += gscore - sc0;
		}
	} 
	a->w = aw[0] > aw[1]? aw[0] : aw[1];
}

void ref_extension(mem_alnreg_t reg, const mem_alnreg_t* reg_dut, const ext_param_t* cur_param, const mem_opt_t* opt) {
	mem_alnreg_t* a = &reg;
	int i;
	int aw[2];
	aw[0] = aw[1] = opt->w;
	int max_off[2];
	if (cur_param->leftQlen > 0) { // left extension
		int qle, tle, gtle, gscore;
		for (i = 0; i < MAX_BAND_TRY; ++i) {
			int prev = a->score;
			aw[0] = opt->w << i;
			a->score = ksw_extend2(cur_param->leftQlen, cur_param->leftQs, cur_param->leftRlen, cur_param->leftRs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, cur_param->h0, &qle, &tle, &gtle, &gscore, &max_off[0]);
			if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
		}
		if (bwa_verbose >= 4) printf("left_ext: qle=%d, tle=%d, gtle=%d, gscore=%d, max_off[0]=%d, score=%d\n", qle, tle, gtle, gscore, max_off[0], a->score);
		if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
			a->qb = cur_param->qBeg - qle, a->rb = cur_param->rbeg_offset - tle;
			a->truesc = a->score;
		} else { // to-end extension
			a->qb = 0, a->rb = cur_param->rbeg_offset - gtle;
			a->truesc = gscore;
		}
	} 

	if (cur_param->rightQlen > 0) { // right extension
		int qle, tle, gtle, gscore, sc0 = a->score;
		for (i = 0; i < MAX_BAND_TRY; ++i) {
			int prev = a->score;
			aw[1] = opt->w << i;
			a->score = ksw_extend2(cur_param->rightQlen, cur_param->rightQs, cur_param->rightRlen, cur_param->rightRs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
			if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
		}
		if (bwa_verbose >= 4) printf("right_ext: qle=%d, tle=%d, gtle=%d, gscore=%d, max_off[1]=%d, score=%d\n", qle, tle, gtle, gscore, max_off[1], a->score);
		// similar to the above
		if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
			a->qe = cur_param->qe_offset + qle, a->re = cur_param->re_offset + tle;
			a->truesc += a->score - sc0;
		} else { // to-end extension
			a->qe = cur_param->qe_offset + cur_param->rightQlen, a->re = cur_param->re_offset + gtle;
			a->truesc += gscore - sc0;
		}
	} 
	a->w = aw[0] > aw[1]? aw[0] : aw[1];
	assert(reg_dut->qb == a->qb);
	assert(reg_dut->qe == a->qe);
	assert(reg_dut->rb == a->rb);
	assert(reg_dut->re == a->re);
	assert(reg_dut->score == a->score);
	assert(reg_dut->truesc == a->truesc);
}

void mem_chain2aln_param(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av, ext_param_t* cur_param, sw_pre_result_t* pre_result, int seed_idx, int batch_idx)
{
	if (bwa_verbose >= 4) printf("[debug] enter mem_chain2aln_param\n");
	int i, k = seed_idx;
	int aw[2]; // aw: actual bandwidth used in extension
	const mem_seed_t *s;
	int64_t tmp;
	uint64_t* srt = pre_result->srt;
	uint8_t* rseq = pre_result->rseq;
	int64_t* rmax = pre_result->rmax;
	cur_param->valid = 0;

	mem_alnreg_t *a;
	s = &c->seeds[(uint32_t)srt[k]];

	if (bwa_verbose >= 4) printf("[debug] start pre-testing\n");
	for (i = 0; i < av->n; ++i) { // test whether extension has been made before
		mem_alnreg_t *p = &av->a[i];
		int64_t rd;
		int qd, w, max_gap;
		if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
		// qd: distance ahead of the seed on query; rd: on reference
		qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
		max_gap = cal_max_gap(opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
		w = max_gap < opt->w? max_gap : opt->w; // bounded by the band width
		if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
		// similar to the previous four lines, but this time we look at the region behind
		qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
		max_gap = cal_max_gap(opt, qd < rd? qd : rd);
		w = max_gap < opt->w? max_gap : opt->w;
		if (qd - rd < w && rd - qd < w) break;
	}
	if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
		if (bwa_verbose >= 4)
			printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment. Confirming whether extension is needed...\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg);
		for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
			const mem_seed_t *t;
			if (srt[i] == 0) continue;
			t = &c->seeds[(uint32_t)srt[i]];
			if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
			if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
			if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
		}
		if (i == c->n) { // no overlapping seeds; then skip extension
			srt[k] = 0; // mark that seed extension has not been performed
			// [QA] Continue to the next seed: continue;
			return;
		}
		if (bwa_verbose >= 4)
			printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
	}
	if (bwa_verbose >= 4) printf("[debug] finish pre-testing\n");

	a = kv_pushp(mem_alnreg_t, *av);
	memset(a, 0, sizeof(mem_alnreg_t));
	a->w = aw[0] = aw[1] = opt->w;
	// [QA] Unnecessary: a->score = a->truesc = -1;
	a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;
	a->qe = l_query, a->re = s->rbeg + s->len;
	
	if (s->qbeg || (s->qbeg + s->len != l_query)) {
		cur_param->valid = 1;
		cur_param->reg = a;
		free(cur_param->leftQs);
		cur_param->leftQs = NULL;
		cur_param->leftQlen = 0;
		cur_param->rightQs = NULL;
		cur_param->rightQlen = 0;
		free(cur_param->leftRs);
		cur_param->leftRs = NULL;
		cur_param->leftRlen = 0;
		cur_param->rightRs = NULL;
		cur_param->rightRlen = 0;
		cur_param->h0 = s->len * opt->a;
		cur_param->regScore = a->score;
		cur_param->qBeg = s->qbeg;
		cur_param->idx = batch_idx; // Batch Index
		cur_param->qe_offset = s->qbeg + s->len;
		cur_param->rbeg_offset = s->rbeg;
		cur_param->re_offset = s->rbeg + s->len;
	} 

	if (bwa_verbose >= 4) printf("[debug] start get the left part\n");
	if (bwa_verbose >= 4) err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] <---\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg);
	if (s->qbeg) { // left extension
		uint8_t *rs, *qs;
		qs = malloc(s->qbeg);
		for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
		tmp = s->rbeg - rmax[0];
		rs = malloc(tmp);
		for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];

		// [QA] Initialize cur_param
		cur_param->leftQlen = s->qbeg;
		cur_param->leftQs = qs;
		cur_param->leftRlen = tmp;
		cur_param->leftRs = rs;
	} 

	if (bwa_verbose >= 4) printf("[debug] start get the right part\n");
	if (s->qbeg + s->len != l_query) { // right extension
		int qe, re;
		qe = s->qbeg + s->len;
		re = s->rbeg + s->len - rmax[0];
		assert(re >= 0);

		cur_param->rightQlen = l_query-qe;
		cur_param->rightQs = query+qe;
		cur_param->rightRlen = rmax[1]-rmax[0]-re;
		cur_param->rightRs = rseq+re;
	} 
	if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);
	cur_param->string_size = ((((cur_param->leftQlen + cur_param->leftRlen 
						    + cur_param->rightQlen + cur_param->rightRlen + 1) >> 1) + 3) >> 2) << 2;

	if (cur_param->valid == 0) {
		if (bwa_verbose >= 4) printf("[debug] optional: calculate seedcov while no parameter");
		// compute seedcov
		for (i = 0, a->seedcov = 0; i < c->n; ++i) {
			const mem_seed_t *t = &c->seeds[i];
			if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
				a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
		}
	}
	if (bwa_verbose >= 4) printf("[debug] return from mem_chain2aln_param to its caller");
}

void mem_chain2aln_seedcov(const mem_chain_t *c, mem_alnreg_t *a)
{
	int i;
	// compute seedcov
	for (i = 0, a->seedcov = 0; i < c->n; ++i) {
		const mem_seed_t *t = &c->seeds[i];
		if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
			a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
	}
	// [QA] Done by Extension: a->w = aw[0] > aw[1]? aw[0] : aw[1];
	// [QA] CRITICAL: free(srt); free(rseq);
}

void mem_chain2aln(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av)
{
	int i, k, max_off[2], aw[2]; // aw: actual bandwidth used in extension
	int64_t rlen, rmax[2], tmp, max = 0;
	const mem_seed_t *s;
	uint8_t *rseq = 0;
	uint64_t *srt;

	if (c->n == 0) return;
	// get the max possible span
	rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) {
		int64_t b, e;
		const mem_seed_t *t = &c->seeds[i];
		b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
		e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
		rmax[0] = rmax[0] < b? rmax[0] : b;
		rmax[1] = rmax[1] > e? rmax[1] : e;
		if (t->len > max) max = t->len;
	}
	rmax[0] = rmax[0] > 0? rmax[0] : 0;
	rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
	if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
		if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
		else rmax[0] = l_pac;
	}
	// retrieve the reference sequence
	rseq = bns_get_seq(l_pac, pac, rmax[0], rmax[1], &rlen);
	assert(rlen == rmax[1] - rmax[0]);

	srt = malloc(c->n * 8);
	for (i = 0; i < c->n; ++i)
		srt[i] = (uint64_t)c->seeds[i].len<<32 | i;
	ks_introsort_64(c->n, srt);

	for (k = c->n - 1; k >= 0; --k) {
		mem_alnreg_t *a;
		s = &c->seeds[(uint32_t)srt[k]];

		for (i = 0; i < av->n; ++i) { // test whether extension has been made before
			mem_alnreg_t *p = &av->a[i];
			int64_t rd;
			int qd, w, max_gap;
			if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
			// qd: distance ahead of the seed on query; rd: on reference
			qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
			max_gap = cal_max_gap(opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
			w = max_gap < opt->w? max_gap : opt->w; // bounded by the band width
			if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
			// similar to the previous four lines, but this time we look at the region behind
			qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
			max_gap = cal_max_gap(opt, qd < rd? qd : rd);
			w = max_gap < opt->w? max_gap : opt->w;
			if (qd - rd < w && rd - qd < w) break;
		}
		if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
			if (bwa_verbose >= 4)
				printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment. Confirming whether extension is needed...\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg);
			for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
				const mem_seed_t *t;
				if (srt[i] == 0) continue;
				t = &c->seeds[(uint32_t)srt[i]];
				if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
				if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
				if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
			}
			if (i == c->n) { // no overlapping seeds; then skip extension
				srt[k] = 0; // mark that seed extension has not been performed
				continue;
			}
			if (bwa_verbose >= 4)
				printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
		}

		a = kv_pushp(mem_alnreg_t, *av);
		memset(a, 0, sizeof(mem_alnreg_t));
		a->w = aw[0] = aw[1] = opt->w;
		a->score = a->truesc = -1;

		if (bwa_verbose >= 4) err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] <---\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg);
		if (s->qbeg) { // left extension
			uint8_t *rs, *qs;
			int qle, tle, gtle, gscore;
			qs = malloc(s->qbeg);
			for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
			tmp = s->rbeg - rmax[0];
			rs = malloc(tmp);
			for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[0] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Left ref:   "); for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
					printf("*** Left query: "); for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
				}
				a->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
				if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
				if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
			}
			// check whether we prefer to reach the end of the query
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
				a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
				a->truesc = a->score;
			} else { // to-end extension
				a->qb = 0, a->rb = s->rbeg - gtle;
				a->truesc = gscore;
			}
			free(qs); free(rs);
		} else a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;

		if (s->qbeg + s->len != l_query) { // right extension
			int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
			qe = s->qbeg + s->len;
			re = s->rbeg + s->len - rmax[0];
			assert(re >= 0);
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[1] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Right ref:   "); for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
					printf("*** Right query: "); for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
				}
				a->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
				if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
				if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
			}
			// similar to the above
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
				a->qe = qe + qle, a->re = rmax[0] + re + tle;
				a->truesc += a->score - sc0;
			} else { // to-end extension
				a->qe = l_query, a->re = rmax[0] + re + gtle;
				a->truesc += gscore - sc0;
			}
		} else a->qe = l_query, a->re = s->rbeg + s->len;
		if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);

		// compute seedcov
		for (i = 0, a->seedcov = 0; i < c->n; ++i) {
			const mem_seed_t *t = &c->seeds[i];
			if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
				a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
		}
		a->w = aw[0] > aw[1]? aw[0] : aw[1];
	}
	free(srt); free(rseq);
}

/*****************************
 * Basic hit->SAM conversion *
 *****************************/

static inline int infer_bw(int l1, int l2, int score, int a, int q, int r)
{
	int w;
	if (l1 == l2 && l1 * a - score < (q + r - a)<<1) return 0; // to get equal alignment length, we need at least two gaps
	w = ((double)((l1 < l2? l1 : l2) * a - score - q) / r + 2.);
	if (w < abs(l1 - l2)) w = abs(l1 - l2);
	return w;
}

static inline int get_rlen(int n_cigar, const uint32_t *cigar)
{
	int k, l;
	for (k = l = 0; k < n_cigar; ++k) {
		int op = cigar[k]&0xf;
		if (op == 0 || op == 2)
			l += cigar[k]>>4;
	}
	return l;
}

void mem_aln2sam(const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_)
{
	int i;
	mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

	if (m_) mtmp = *m_, m = &mtmp;
	// set flag
	p->flag |= m? 0x1 : 0; // is paired in sequencing
	p->flag |= p->rid < 0? 0x4 : 0; // is mapped
	p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
	if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
		p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
	if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
		m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
	p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
	p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand

	// print up to CIGAR
	kputs(s->name, str); kputc('\t', str); // QNAME
	kputw((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc('\t', str); // FLAG
	if (p->rid >= 0) { // with coordinate
		kputs(bns->anns[p->rid].name, str); kputc('\t', str); // RNAME
		kputl(p->pos + 1, str); kputc('\t', str); // POS
		kputw(p->mapq, str); kputc('\t', str); // MAPQ
		if (p->n_cigar) { // aligned
			for (i = 0; i < p->n_cigar; ++i) {
				int c = p->cigar[i]&0xf;
				if (c == 3 || c == 4) c = which? 4 : 3; // use hard clipping for supplementary alignments
				kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
			}
		} else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
	} else kputsn("*\t0\t0\t*", 7, str); // without coordinte
	kputc('\t', str);

	// print the mate position if applicable
	if (m && m->rid >= 0) {
		if (p->rid == m->rid) kputc('=', str);
		else kputs(bns->anns[m->rid].name, str);
		kputc('\t', str);
		kputl(m->pos + 1, str); kputc('\t', str);
		if (p->rid == m->rid) {
			int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
			int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
			if (m->n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
			else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
		} else kputc('0', str);
	} else kputsn("*\t0\t0", 5, str);
	kputc('\t', str);

	// print SEQ and QUAL
	if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
		kputsn("*\t*", 3, str);
	} else if (!p->is_rev) { // the forward strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar) {
			if (which && ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3)) qb += p->cigar[0]>>4;
			if (which && ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3)) qe -= p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	} else { // the reverse strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar) {
			if (which && ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3)) qe -= p->cigar[0]>>4;
			if (which && ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3)) qb += p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	}

	// print optional tags
	if (p->n_cigar) {
		kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
		kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
	}
	if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
	if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(p->sub, str); }
	if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
	if (!(p->flag & 0x100)) { // not multi-hit
		for (i = 0; i < n; ++i)
			if (i != which && !(list[i].flag&0x100)) break;
		if (i < n) { // there are other primary hits; output them
			kputsn("\tSA:Z:", 6, str);
			for (i = 0; i < n; ++i) {
				const mem_aln_t *r = &list[i];
				int k;
				if (i == which || (list[i].flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hit
				kputs(bns->anns[r->rid].name, str); kputc(',', str);
				kputl(r->pos+1, str); kputc(',', str);
				kputc("+-"[r->is_rev], str); kputc(',', str);
				for (k = 0; k < r->n_cigar; ++k) {
					kputw(r->cigar[k]>>4, str); kputc("MIDSH"[r->cigar[k]&0xf], str);
				}
				kputc(',', str); kputw(r->mapq, str);
				kputc(',', str); kputw(r->NM, str);
				kputc(';', str);
			}
		}
	}
	if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
	kputc('\n', str);
}

/************************
 * Integrated interface *
 ************************/

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a)
{
	int mapq, l, sub = a->sub? a->sub : opt->min_seed_len * opt->a;
	double identity;
	sub = a->csub > sub? a->csub : sub;
	if (sub >= a->score) return 0;
	l = a->qe - a->qb > a->re - a->rb? a->qe - a->qb : a->re - a->rb;
	identity = 1. - (double)(l * opt->a - a->score) / (opt->a + opt->b) / l;
	if (a->score == 0) {
		mapq = 0;
	} else if (opt->mapQ_coef_len > 0) {
		double tmp;
		tmp = l < opt->mapQ_coef_len? 1. : opt->mapQ_coef_fac / log(l);
		tmp *= identity * identity;
		mapq = (int)(6.02 * (a->score - sub) / opt->a * tmp * tmp + .499);
	} else {
		mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499);
		mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
	}
	if (a->sub_n > 0) mapq -= (int)(4.343 * log(a->sub_n+1) + .499);
	if (mapq > 60) mapq = 60;
	if (mapq < 0) mapq = 0;
	return mapq;
}

// TODO (future plan): group hits into a uint64_t[] array. This will be cleaner and more flexible
void mem_reg2sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m)
{
	kstring_t str;
	kvec_t(mem_aln_t) aa;
	int k;

	kv_init(aa);
	str.l = str.m = 0; str.s = 0;
	for (k = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		mem_aln_t *q;
		if (p->score < opt->T) continue;
		if (p->secondary >= 0 && !(opt->flag&MEM_F_ALL)) continue;
		if (p->secondary >= 0 && p->score < a->a[p->secondary].score * .5) continue;
		q = kv_pushp(mem_aln_t, aa);
		*q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
		q->flag |= extra_flag; // flag secondary
		if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
		if (k && p->secondary < 0) // if supplementary
			q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
		if (k && q->mapq > aa.a[0].mapq) q->mapq = aa.a[0].mapq;
	}
	if (aa.n == 0) { // no alignments good enough; then write an unaligned record
		mem_aln_t t;
		t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
		t.flag |= extra_flag;
		mem_aln2sam(bns, &str, s, 1, &t, 0, m);
	} else {
		for (k = 0; k < aa.n; ++k)
			mem_aln2sam(bns, &str, s, aa.n, aa.a, k, m);
		for (k = 0; k < aa.n; ++k) free(aa.a[k].cigar);
		free(aa.a);
	}
	s->sam = str.s;
}

/*
void retrieve_output_memory(const ext_param_t* param_batch, const int left_idx, const int batch_idx, const int8_t* output_space, const mem_opt_t* opt) {
	int i;
	int16_t* p_res = (int16_t*)output_space;
	for (i=left_idx; i<=batch_idx; i++) {
		if (param_batch[i].valid) {
			mem_alnreg_t ref_reg = *(param_batch[i].reg);
			param_batch[i].reg->qb = *p_res;
			p_res++;
			param_batch[i].reg->qe = *p_res + param_batch[i].qe_offset;
			p_res++;
			param_batch[i].reg->rb = *p_res + param_batch[i].rbeg_offset;
			p_res++;
			param_batch[i].reg->re = *p_res + param_batch[i].re_offset;
			p_res++;
			param_batch[i].reg->score = *p_res;
			p_res++;
			param_batch[i].reg->truesc = *p_res;
			p_res++;
			param_batch[i].reg->w = *p_res;
			p_res++;
			p_res++;
			if (bwa_verbose >= 4) ref_extension(ref_reg, param_batch[i].reg, &param_batch[i], opt);
		}
	}
}*/

void retrieve_output_memory(const ext_param_t* param_batch, const int8_t* output_space, const int taskNum, int start) {
	printf("\n[retrieve_output_memory] #task = %d, start = %d\n", taskNum, start);
        int i;
        int16_t* p_res = (int16_t*)output_space;
	int32_t cur_idx = -1;
        for (i=0; i<taskNum; i++) {
		cur_idx = *((int32_t*)(p_res));
		cur_idx -= start;
		printf ("[retrieve_output_memory] the index of the %d-th task is: %d\n", i, cur_idx);
		assert (param_batch[cur_idx].valid);
		p_res += 2;
                param_batch[cur_idx].reg->qb = *p_res;
                p_res++;
                param_batch[cur_idx].reg->qe = *p_res + param_batch[cur_idx].qe_offset;
                p_res++;
                param_batch[cur_idx].reg->rb = *p_res + param_batch[cur_idx].rbeg_offset;
                p_res++;
                param_batch[cur_idx].reg->re = *p_res + param_batch[cur_idx].re_offset;
                p_res++;
                param_batch[cur_idx].reg->score = *p_res;
                p_res++;
                param_batch[cur_idx].reg->truesc = *p_res;
                p_res++;
                param_batch[cur_idx].reg->w = *p_res;
                p_res++;
                p_res++;
        }
		printf ("\n[retrieve_output_memory] finish retrieving\n");
}

/*
void retrieve_output_memory(const ext_param_t* param_batch, const int16_t* output_space, const int taskNum, int start) {
	printf("\n[retrieve_output_memory] #task = %d\n", taskNum);

	int i, j;

	int32_t cur_idx = -1;
	int32_t  valid_cache_line;
	int16_t* p_res;
	int16_t* p_base;

	valid_cache_line = taskNum * 5 / 16 + 1;
	p_res = (int16_t*)malloc( sizeof(int16_t) * valid_cache_line * 32 );
	p_base = (int16_t*) output_space;

	for (j = 0; j < valid_cache_line; ++j, p_base += 32 ) {
		for (i = 0; i < 32; i++) {
			p_res[j * 32 + i] = *(int16_t *)(p_base + (31 - i));
		}
	}

	// int16_t* p_res = (int16_t*)output_space;
	for (i = 0; i < taskNum; i++) {
		// cur_idx = *((int32_t*)(p_res));
		cur_idx = (int32_t)(p_res[i * 10 + 1]);
		cur_idx -= start;
		// p_res += 3;
		param_batch[cur_idx].reg->qb = p_res[i * 10 + 3];
		// p_res--;
		param_batch[cur_idx].reg->qe = p_res[i * 10 + 2] + param_batch[cur_idx].qe_offset;
		// p_res += 3;
		param_batch[cur_idx].reg->rb = p_res[i * 10 + 5] + param_batch[cur_idx].rbeg_offset;
		// p_res--;
		param_batch[cur_idx].reg->re = p_res[i * 10 + 4] + param_batch[cur_idx].re_offset;
		// p_res += 3;
		param_batch[cur_idx].reg->score = p_res[i * 10 + 7];
		// p_res--;
		param_batch[cur_idx].reg->truesc = p_res[i * 10 + 6];
		// p_res += 3;
		param_batch[cur_idx].reg->w = p_res[i * 10 + 9];
		// p_res++;
	}

	free(p_res);

	printf("done\n");
}
*/

void fill_input_memory(const ext_param_t* param_batch, const int left_idx, const int batch_idx, int8_t* input_space, const mem_opt_t* opt) {
	int i,j;
	int param_idx = 0;
	int numOfReads = *((int*)(input_space+8));
	printf("\n[fill_input_memory] #task = %d\n", numOfReads);
	int8_t* p_param = input_space + 32; 
	int8_t* p_string = input_space + 32 + 32 * numOfReads;
	int32_t offset = 8 + numOfReads * 8;

	int16_t leftMaxIns, leftMaxDel, rightMaxIns, rightMaxDel;

	for (i=left_idx; i<=batch_idx; i++) {
		if (param_batch[i].valid) {
			*((int16_t*)(&p_param[0])) = param_batch[i].leftQlen;
			*((int16_t*)(&p_param[2])) = param_batch[i].leftRlen;
			*((int16_t*)(&p_param[4])) = param_batch[i].rightQlen;
			*((int16_t*)(&p_param[6])) = param_batch[i].rightRlen;
			*((int32_t*)(&p_param[8])) = offset;
			offset += (param_batch[i].string_size >> 2);
			*((int16_t*)(&p_param[12])) = param_batch[i].regScore;
			*((int16_t*)(&p_param[14])) = param_batch[i].qBeg;
			*((int16_t*)(&p_param[16])) = param_batch[i].h0;
			// [QA] no p_param[18...19]
			// [QA] These four variables only fit for the default matrix
			leftMaxIns = (param_batch[i].leftQlen + opt->pen_clip5 - opt->o_ins) / opt->e_ins + 1;
			leftMaxDel = (param_batch[i].leftQlen + opt->pen_clip5 - opt->o_del) / opt->e_del + 1;
			rightMaxIns = (param_batch[i].rightQlen + opt->pen_clip3 - opt->o_ins) / opt->e_ins + 1;
			rightMaxDel = (param_batch[i].rightQlen + opt->pen_clip3 - opt->o_del) / opt->e_del + 1;
			*((int16_t*)(&p_param[20])) = leftMaxIns;
			*((int16_t*)(&p_param[22])) = leftMaxDel;
			*((int16_t*)(&p_param[24])) = rightMaxIns;
			*((int16_t*)(&p_param[26])) = rightMaxDel;
			*((int32_t*)(&p_param[28])) = param_batch[i].idx;
			printf("\n[fill_input_memory] the %d-th param has index: %d\n", i, param_batch[i].idx);

			int bp_num = 0;
			uint32_t* bp_val = (uint32_t*)p_string;
			for (j=0; j<param_batch[i].leftQlen; j++) {
				bp_num ++;
				*bp_val = ((*bp_val) << 4) | param_batch[i].leftQs[j];
				if (bp_num % 8 == 0) bp_val = bp_val + 1;
			}
			for (j=0; j<param_batch[i].rightQlen; j++) {
				bp_num ++;
				*bp_val = ((*bp_val) << 4) | param_batch[i].rightQs[j];
				if (bp_num % 8 == 0) bp_val = bp_val + 1;
			}
			for (j=0; j<param_batch[i].leftRlen; j++) {
				bp_num ++;
				*bp_val = ((*bp_val) << 4) | param_batch[i].leftRs[j];
				if (bp_num % 8 == 0) bp_val = bp_val + 1;
			}
			for (j=0; j<param_batch[i].rightRlen; j++) {
				bp_num ++;
				*bp_val = ((*bp_val) << 4) | param_batch[i].rightRs[j];
				if (bp_num % 8 == 0) bp_val = bp_val + 1;
			}
			assert (bp_num == param_batch[i].leftQlen + param_batch[i].leftRlen + param_batch[i].rightQlen + param_batch[i].rightRlen);
			while (bp_num % 8) *bp_val = ((*bp_val) << 4), bp_num++;
			p_param += 32;
			p_string += param_batch[i].string_size;
			param_idx++;
		}
	}
	assert(p_param == input_space + 32 + 32 * numOfReads);
	assert(numOfReads == param_idx);
}

void sw_extend(unsigned short qs_baddr, int *qs, unsigned short ts_baddr, short qlen, short tlen, char o_ins,
			   char e_ins, char o_del, char e_del, char penClip, char w_in, char h0, short *regScore, short qBeg, short max_ins, short max_del,
			   short *w_ret, short *qle_ret, short *tle_ret, short *gtle_ret, short *gscore_ret, short *maxoff_ret)
{
	//#pragma HLS INLINE

	int i, j;
	int k, l;
	short max_i, max_ie, max_off;
	short gscore;
	char max_j;
	char oe_del = o_del + e_del;
	char oe_ins = o_ins + e_ins;
	short beg, end;
	char backw_tmp=0;
	char backw_reg=0;
	char forw_update=0;
	char forw_tmp=0;
	char forw_reg=0;
	short abs_mj_m_i;
	char tmp_ehh_m_eins;
	char tmp_eme;
	char h1_init_val;
	char max;
	char h, e;
	char e_tmp;
	char h_tmp;
	char h1_reg;
	char t, f = 0, h1, m = 0;
	char mj = -1;
	char q_i = 0, q_j = 0;
	short prev;
	char isBreak;
	char aw1;
	char aw_tmp;
	char h0_arr[2];
#pragma HLS ARRAY_PARTITION variable=h0_arr complete dim=0
	//	short sc0;
	//	short h0;

	const char my_mat[5][5]={{1, -4, -4, -4, -1}, {-4, 1, -4, -4, -1}, {-4, -4, 1, -4, -1}, {-4, -4, -4, 1, -1}, {-1, -1, -1, -1, -1}};
#pragma HLS ARRAY_PARTITION variable=my_mat complete dim=0
	char eh_h [256];
#pragma HLS ARRAY_MAP variable=eh_h instance=eh_arr vertical
#pragma HLS RESOURCE variable=eh_h core=RAM_2P_BRAM
	char eh_e [256];
#pragma HLS ARRAY_MAP variable=eh_e instance=eh_arr vertical
#pragma HLS RESOURCE variable=eh_e core=RAM_2P_BRAM

	//	sc0 = *regScore;
	//	h0_arr[0] = h0_ori;
	//	h0_arr[1] = sc0;
	//	qle = -1;
	//	tle = -1;
	//	gtle = -1;
	max = h0;
	max_i = max_j = -1;
	max_ie = -1;
	gscore = -1;
	max_off = 0;

	k = 0;
	isBreak = 0;
ext_while_loop : while ((k < 2) && (!isBreak))
				 {
#pragma HLS LOOP_TRIPCOUNT min=2 max=2
					 prev = *regScore;
					 aw_tmp = w_in << k;
					 aw1 = aw_tmp < max_ins ? aw_tmp : max_ins;
					 aw1 = aw1 < max_del ? aw1 : max_del;
					 beg = 0;
					 end = qlen;
					 if (h0 > oe_ins) {
						 tmp_eme = h0 - oe_ins;
					 }
					 else {
						 tmp_eme = 0;
					 }
					 h1_init_val = h0 - o_del;
target_loop : for (i = 0; i < tlen; i++) {
					 f = 0; m = 0; mj = -1;
					 //while(pe_seeds.empty());
					 //			q_i = ts[ts_baddr + i];
					 q_i = qs[ts_baddr + i];
					 h1_init_val -= e_del;
					 h1 = h1_init_val;
					 if (h1 < 0) h1 = 0;
					 if (beg < i - aw1) beg = i - aw1;
					 if (end > i + aw1 + 1) end = i + aw1 + 1;
					 if (end > qlen) end = qlen;
					 backw_tmp = 0; backw_reg = 0;
					 forw_tmp = 0; forw_reg = 0;
					 forw_update = 0;
query_loop : for (j = beg; j < end; ++j) {
#pragma HLS LOOP_TRIPCOUNT min=50 max=50
#pragma AP pipeline II=1
#pragma AP dependence variable=eh_e array inter false
#pragma AP dependence variable=eh_h array inter false
					 q_j = qs[qs_baddr + j];
					 h_tmp = eh_h[j];// get H(i-1,j-1) and E(i-1,j)
					 e_tmp = eh_e[j];
					 if (i == 0) {
						 e = 0;
						 if (j == 0) {
							 h = h0;
						 }
						 else if (j == 1) {
							 h = tmp_eme;
						 }
						 else {
							 tmp_eme -= e_ins;
							 if (tmp_eme > 0) {
								 h = tmp_eme;
							 }
							 else {
								 h = 0;
							 }
						 }
					 }
					 else {
						 e = e_tmp;
						 h = h_tmp;
					 }
					 h1_reg = h1;
					 h += my_mat[q_i][q_j];
					 h = h > e? h : e;
					 h = h > f? h : f;
					 h1 = h;             // save H(i,j) to h1 for the next column
					 if (h1_reg == 0) {
						 backw_tmp = 0;
					 }
					 else {
						 backw_tmp++;
					 }
					 if (m <= h)
					 {
						 mj = j;
						 m = h;
						 backw_reg = backw_tmp;
					 }
					 if (j >= mj+2) {
						 if (forw_update == 0) { //((h1_reg == 0) &&
							 if (h1_reg == 0) {
								 forw_update = 1;
							 }
							 else {
								 forw_tmp++;
							 }
						 }
					 }
					 else {
						 forw_tmp = 0;
						 forw_update = 0;
					 }
					 t = h - oe_del;
					 t = t > 0? t : 0;
					 e -= e_del;
					 e = e > t? e : t;   // computed E(i+1,j)
					 eh_e[j] = e; // save E(i+1,j) for the next row
					 eh_h[j] = h1_reg;          // set H(i,j-1) for the next row
					 t = h - oe_ins;
					 t = t > 0? t : 0;
					 f -= e_ins;
					 f = f > t? f : t;   // computed F(i,j+1)
			 }
			 eh_h[end] = h1;
			 eh_e[end] = 0;
			 if ((forw_update == 0) && (h1 != 0)) {
				 if ((j >= mj+2) || (forw_tmp != 0)) {
					 forw_tmp++;
				 }
			 }
			 //			fprintf(fp_test, "\t k = %d;\t i = %d;\t j = %d;\t beg = %d;\t end = %d;\t gscore = %d;\t h1 = %d;\t max_ie = %d;\t", k, i, j, beg, end, gscore, h1, max_ie);
			 if (j == qlen) {
				 max_ie = gscore > h1? max_ie : i;
				 gscore = gscore > h1? gscore : h1;
			 }
			 //			fprintf(fp_test, "\t m = %d;\t max = %d;\t mj = %d;\t max_j = %d;\t\n ", m, max, mj, max_j);
			 if (m == 0) break;
			 if (m > max) {
				 max = m; max_i = i; max_j = mj;
				 if (mj >= i) abs_mj_m_i = mj - i;
				 else abs_mj_m_i = i - mj;
				 if (max_off < abs_mj_m_i) max_off = abs_mj_m_i;
			 }
			 // update beg and end for the next round
			 //for (j = mj; j >= beg && eh_h[j]; --j);
			 //beg = j + 1;
			 //for (j = mj + 2; j <= end && eh_h[j]; ++j);
			 //end = j;
			 j = mj - backw_reg;
			 beg = j + 1;
			 j = mj + 2 + forw_tmp;
			 end = j;
			  }
			  *qle_ret = max_j + 1;
			  *tle_ret = max_i + 1;
			  *gtle_ret = max_ie + 1;
			  *gscore_ret = gscore;
			  *maxoff_ret = max_off;
			  *regScore = max;
			  if (max == prev || ( max_off < (aw_tmp >> 1) + (aw_tmp >> 2))) isBreak = 1;
			  k++;
				 }
				 *w_ret = aw_tmp;
				 //	*score = max;
}

void extension(int8_t* input_space, int8_t* output_space)
{   
	int * a; //cmost_malloc_1d( &a, "total_input.dat" , 4, 409600);
	a = (int*)input_space;

    	int * c; //cmost_malloc_1d( &c, "zero.dat" , 4, 8192);
	c = (int*)output_space;

#pragma cmost graph_begin
#pragma cmost task_begin name="vec_add" lib="true"
		int ii = 0;
		int idx = 0;
		char o_del;
		char e_del;
		char o_ins;
		char e_ins;
		char penClip[2];
#pragma HLS ARRAY_PARTITION variable=penClip complete dim=0
		char w_in;
		int taskNums;
		int tmp_compar;
		tmp_compar = a[0];
		o_del = tmp_compar & 0xFF; //6
		e_del = (tmp_compar >> 8) & 0xFF; //1
		o_ins = (tmp_compar >> 16) & 0xFF; //6
		e_ins = (tmp_compar >> 24) & 0xFF; //1
		tmp_compar = a[1];
		penClip[0] = tmp_compar & 0xFF; //100
		penClip[1] = (tmp_compar >> 8) & 0xFF;
		w_in = (tmp_compar >> 16) & 0xFF;
		//h0 = (tmp_compar >> 24) & 0xFF; //45
		taskNums = a[2];

		for (ii=0; ii<taskNums; ii++)
		{
			int seed_index;
			short qlen[2];
#pragma HLS ARRAY_PARTITION variable=qlen complete dim=0
			short tlen[2];
#pragma HLS ARRAY_PARTITION variable=tlen complete dim=0
			short max_ins[2];
#pragma HLS ARRAY_PARTITION variable=max_ins complete dim=0
			short max_del[2];
#pragma HLS ARRAY_PARTITION variable=max_del complete dim=0
			char h0;
			short regScore;
			short qBeg_ori;

			short qle;
			short tle;
			short gtle;
			short gscore;
			short maxoff;
			short qBeg;
			short rBeg;
			short qEnd;
			short rEnd;
			short score;
			short trueScore;
			short width;

			unsigned char i;
			short k, l;
			short qlen2;
			short tlen2;
			unsigned short qs_baddr;
			unsigned short ts_baddr;
			short aw[2];
#pragma HLS ARRAY_PARTITION variable=aw complete dim=0
			int tmp_parame;
			int qr_offset;
			short sc0;
			short h0_arr[2];
#pragma HLS ARRAY_PARTITION variable=h0_arr complete dim=0

			int qrLen_div8;
			int tmp_qr_data;
			int tmp_query_mem[256];
#pragma HLS RESOURCE variable=tmp_query_mem core=RAM_2P_BRAM
			int query_mem[2048];
#pragma HLS RESOURCE variable=query_mem core=RAM_2P_BRAM
			//		int target_mem[2048];
			//	#pragma HLS RESOURCE variable=target_mem core=RAM_2P_BRAM
			int param_mem[8];

			//memcpy(param_mem, &a[(ii+1)*8], (8)*4);
			for (idx=0; idx<8; idx++)
			{
				param_mem[idx] = a[(ii+1)*8 + idx];
			}
			//tmp_parame = totalinp[ii*8 + 0];
			//while(pe_seeds.empty());
			//pe_seeds.read(seed_index);
			//tmp_parame = totalinp[ii*8 + 1];
			//while(pe_seeds.empty());
			//pe_seeds.read(tmp_parame);
			//		tmp_parame = totalinp[(ii+1)*8 + 0];
			tmp_parame = param_mem[0];
			qlen[0] = tmp_parame & 0xFFFF; //55
			tlen[0] = (tmp_parame >> 16) & 0xFFFF; //105
			//tmp_parame = totalinp[ii*8 + 2];
			//while(pe_seeds.empty());
			//pe_seeds.read(tmp_parame);
			//		tmp_parame = totalinp[(ii+1)*8 + 1];
			tmp_parame = param_mem[1];
			qlen[1] = tmp_parame & 0xFFFF; //55
			tlen[1] = (tmp_parame >> 16) & 0xFFFF; //105
			//		tmp_parame = totalinp[(ii+1)*8 + 2];
			tmp_parame = param_mem[2];
			qr_offset = tmp_parame;
			//		tmp_parame = totalinp[(ii+1)*8 + 3];
			tmp_parame = param_mem[3];
			regScore = tmp_parame & 0xFFFF; //100
			qBeg_ori = (tmp_parame >> 16) & 0xFFFF;
			//		tmp_parame = totalinp[(ii+1)*8 + 4];
			tmp_parame = param_mem[4];
			h0 = tmp_parame & 0xFFFF;
			//max_del[0] = (tmp_parame >> 16) & 0xFFFF;
			//tmp_parame = totalinp[ii*8 + 3];
			//while(pe_seeds.empty());
			//pe_seeds.read(tmp_parame);
			//		tmp_parame = totalinp[ii*8 + 2];
			//		o_del = tmp_parame & 0xFF; //6
			//		e_del = (tmp_parame >> 8) & 0xFF; //1
			//		o_ins = (tmp_parame >> 16) & 0xFF; //6
			//		e_ins = (tmp_parame >> 24) & 0xFF; //1
			//tmp_parame = totalinp[ii*8 + 4];
			//while(pe_seeds.empty());
			//pe_seeds.read(tmp_parame);
			//		tmp_parame = totalinp[ii*8 + 3];
			//		penClip[0] = tmp_parame & 0xFF; //100
			//		penClip[1] = (tmp_parame >> 8) & 0xFF;
			//		w_in = (tmp_parame >> 16) & 0xFF;
			//		h0 = (tmp_parame >> 24) & 0xFF; //45
			//tmp_parame = totalinp[ii*8 + 5];
			//pe_seeds.read(tmp_parame);
			//		tmp_parame = totalinp[(ii+1)*8 + 5];
			tmp_parame = param_mem[5];
			max_ins[0] = tmp_parame & 0xFFFF;
			max_del[0] = (tmp_parame >> 16) & 0xFFFF;
			//		tmp_parame = totalinp[(ii+1)*8 + 6];
			tmp_parame = param_mem[6];
			max_ins[1] = tmp_parame & 0xFFFF;
			max_del[1] = (tmp_parame >> 16) & 0xFFFF;
			//tmp_parame = totalinp[ii*8 + 6];
			//while(pe_seeds.empty());
			//pe_seeds.read(tmp_parame);
			//		tmp_parame = totalinp[(ii+1)*8 + 6];
			//		regScore = tmp_parame & 0xFFFF; //100
			//		qBeg_ori = (tmp_parame >> 16) & 0xFFFF;
			//tmp_parame = totalinp[ii*8 + 7];
			//while(pe_seeds.empty());
			//pe_seeds.read(qr_offset);
			//		tmp_parame = totalinp[(ii+1)*8 + 7];
			//		qr_offset = tmp_parame;

			aw[0] = w_in;
			aw[1] = w_in;
			qBeg = 0;
			qEnd = qlen[1];
			rBeg = 0;
			rEnd = 0;
			trueScore = regScore;
			qle = -1;
			tle = -1;
			gtle = -1;
			gscore = -1;
			maxoff = -1;
			qlen2 = qlen[0] + qlen[1];
			tlen2 = tlen[0] + tlen[1];
			qrLen_div8 = qlen2+tlen2;
			if ((qrLen_div8 & 0x00000007) != 0) {
				qrLen_div8 = (qrLen_div8 >> 3) + 1;
			}
			else {
				qrLen_div8 = (qrLen_div8 >> 3);
			}

			//memcpy(tmp_query_mem, &a[qr_offset], (qrLen_div8)*4);
			for (idx=0; idx<qrLen_div8; idx++)
			{
				tmp_query_mem[idx] = a[qr_offset + idx];
			}

			for (k=0; k<qrLen_div8; k++) {
				tmp_qr_data = tmp_query_mem[k];
				for (l=0; l<8; l++) {
					query_mem[k*8 + l] = (tmp_qr_data & 0xF0000000) >> 28;
					tmp_qr_data <<= 4;
				}
			}
			//	load_query : for (k = 0; k < qlen2; k++) {
			//#pragma HLS PIPELINE II=1
			//#pragma HLS LOOP_TRIPCOUNT min=50 max=50
			//			query_mem[k] = totalinp[qr_offset];
			//			qr_offset++;
			//while(pe_seeds.empty());
			//pe_seeds.read(query_data);
			//query_mem[k] = query_data;
			//		}
			//	memcpy(target_mem, &totalinp[qr_offset+qlen2], tlen2*4);
			//	load_target : for (k = 0; k < tlen2; k++) {
			//#pragma HLS PIPELINE II=1
			//#pragma HLS LOOP_TRIPCOUNT min=100 max=100
			//		    target_mem[k] = totalinp[qr_offset];
			//			qr_offset++;
			//while(pe_seeds.empty());
			//pe_seeds.read(target_data);
			//target_mem[k] = target_data;
			//		}
			qs_baddr = 0;
			//		ts_baddr = 0;
			ts_baddr = qlen2;
left_right_loop : for (i=0; i<2; i++)
				  {
					  //#pragma HLS PIPELINE
#pragma HLS LOOP_TRIPCOUNT min=2 max=2
					  sc0 = regScore;
					  h0_arr[0] = h0;
					  h0_arr[1] = sc0;
					  if (qlen[i] > 0) {
						  sw_extend(qs_baddr, query_mem, ts_baddr, qlen[i], tlen[i], o_ins, e_ins, o_del, e_del, penClip[i],
							  w_in, h0_arr[i], &regScore, qBeg_ori, max_ins[i], max_del[i], &aw[i], &qle, &tle, &gtle, &gscore, &maxoff);
						  score = regScore;
						  if (gscore <= 0 || gscore <= (regScore - penClip[i])) {
							  if (i == 0) {
								  qBeg = qBeg_ori - qle;
								  rBeg = -tle;
								  trueScore = regScore;
							  }
							  else {
								  qEnd = qle;
								  rEnd = tle;
								  trueScore += regScore - sc0;
							  }
						  }
						  else {
							  if (i == 0) {
								  qBeg = 0;
								  rBeg = -gtle;
								  trueScore = gscore;
							  }
							  else {
								  qEnd = qlen[1];
								  rEnd = gtle;
								  trueScore += gscore - sc0;
							  }
						  }
					  }
					  qs_baddr += qlen[i];
					  ts_baddr += tlen[i];
				  }
				  if (aw[0] > aw[1]) width = aw[0];
				  else width = aw[1];

				  //		fprintf(fp_test, "\t ii = %d;\t tle = %d;\t gtle = %d;\t rBeg = %d;\n ", ii, tle, gtle, rBeg);
				  //		fprintf(fp_test, "\t ii = %d;\t qBeg_ori = %d;\t qle = %d;\t qBeg = %d;\n ", ii, qBeg_ori, qle, qBeg);
				  //while(pe_matchs.full());
				  //pe_matchs.write(seed_index);
				  c[ii*4 + 0] = (qBeg & 0xFFFF) | ((qEnd<<16) & 0xFFFF0000);
				  //while(pe_matchs.full());
				  //pe_matchs.write(result_data);
				  c[ii*4 + 1] = (rBeg & 0xFFFF) | ((rEnd<<16) & 0xFFFF0000);
				  //while(pe_matchs.full());
				  //pe_matchs.write(result_data);
				  c[ii*4 + 2] = (score & 0xFFFF) | ((trueScore<<16) & 0xFFFF0000);
				  //while(pe_matchs.full());
				  //pe_matchs.write(result_data);
				  c[ii*4 + 3] = width & 0xFFFF;
		}

#pragma cmost task_end
#pragma cmost graph_end
}

typedef struct {
	int chn_idx;
	int seed_idx;
} coord_t;

// [QA] The batched-processing version of mem_align1_core
mem_alnreg_v* mem_align1_core_batched(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, bseq1_t* seqs, int start, int batch_size)
{
	int i,batch_idx;
	mem_chain_v* chn_batch = (mem_chain_v*)malloc(sizeof(mem_chain_v)*batch_size);
	mem_alnreg_v* regs_batch = (mem_alnreg_v*)malloc(sizeof(mem_alnreg_v)*batch_size);

	for (batch_idx=start; batch_idx<start+batch_size; batch_idx++) {
		for (i = 0; i < seqs[batch_idx].l_seq; ++i) // convert to 2-bit encoding if we have not done so
			seqs[batch_idx].seq[i] = seqs[batch_idx].seq[i] < 4? seqs[batch_idx].seq[i] : nst_nt4_table[(int)(seqs[batch_idx].seq[i])];

		chn_batch[batch_idx-start] = mem_chain(opt, bwt, bns->l_pac, seqs[batch_idx].l_seq, (uint8_t*)(seqs[batch_idx].seq));
		chn_batch[batch_idx-start].n = mem_chain_flt(opt, chn_batch[batch_idx-start].n, chn_batch[batch_idx-start].a);
		if (bwa_verbose >= 4) mem_print_chain(bns, &chn_batch[batch_idx-start]);

		kv_init(regs_batch[batch_idx-start]);
	}

	coord_t* coord_batch = (coord_t*)malloc(sizeof(coord_t)*batch_size);
	sw_pre_result_t* result_batch = (sw_pre_result_t*)malloc(sizeof(sw_pre_result_t)*batch_size);
	memset(result_batch, 0, sizeof(sw_pre_result_t)*batch_size);
	ext_param_t* param_batch = (ext_param_t*)malloc(sizeof(ext_param_t)*batch_size);
	memset(param_batch, 0, sizeof(ext_param_t)*batch_size);

	/* backup	
	// [QA] Input Memory Space
	int8_t input_space[INPUT_SPACE_SIZE];
	int8_t output_space[OUTPUT_SPACE_SIZE];
	input_space[0] = opt->o_del;
	input_space[1] = opt->e_del;
	input_space[2] = opt->o_ins;
	input_space[3] = opt->e_ins;
	input_space[4] = opt->pen_clip5;
	input_space[5] = opt->pen_clip3;
	input_space[6] = opt->w;
	// [QA] input_space[7] = 0;
	int32_t* p_readNo = (int32_t*)(input_space+8);
	// [QA] input_space[12...31]
	*/

	int32_t taskNum;
	int32_t* p_readNo = &taskNum;
	

	// [QA] Initialize the coordinates for all the reads
	// [QA] Chain Idx = 0, Seed Idx = #Seeds - 1
	// [QA] If a read has no chain, then invalidate its parameter structure
	// [QA] If a read has at least one chain, then prepare its first preprocessing result and parameter structures
	for (batch_idx=start; batch_idx<start+batch_size; batch_idx++) {
		coord_batch[batch_idx-start].chn_idx = 0;
		coord_batch[batch_idx-start].seed_idx = chn_batch[batch_idx-start].n>0? chn_batch[batch_idx-start].a[coord_batch[batch_idx-start].chn_idx].n-1:  -1;
	}

	int left_idx = start, right_idx = start+batch_size;
	int filled_space = 32; // 32 bytes for the header
	*p_readNo = 0;
	int start_read_idx = left_idx;
	// [QA] the batch bound shrinks gradually
	while (left_idx < right_idx) {
		if (bwa_verbose >= 4) printf("[debug] start processing a line of seeds\n");
		filled_space = 32; // 32 bytes for the header
		*p_readNo = 0;
		start_read_idx = left_idx;
		// [QA] prepare smith-waterman inputs
		if (bwa_verbose >= 4) printf("[debug] 1st: prepocessing and param retrieving\n");
		for (batch_idx=left_idx; batch_idx<right_idx; batch_idx++) {
			if (chn_batch[batch_idx-start].n > coord_batch[batch_idx-start].chn_idx) {
				if (coord_batch[batch_idx-start].seed_idx >= 0) {
					if (coord_batch[batch_idx-start].seed_idx == chn_batch[batch_idx-start].a[coord_batch[batch_idx-start].chn_idx].n-1) {
						mem_chain2aln_preprocess(opt, bns->l_pac, pac, seqs[batch_idx].l_seq, (uint8_t*)(seqs[batch_idx].seq), &chn_batch[batch_idx-start].a[coord_batch[batch_idx-start].chn_idx], &result_batch[batch_idx-start]);
					}
					mem_chain2aln_param(opt, bns->l_pac, pac, seqs[batch_idx].l_seq, (uint8_t*)(seqs[batch_idx].seq), &chn_batch[batch_idx-start].a[coord_batch[batch_idx-start].chn_idx], &regs_batch[batch_idx-start], &param_batch[batch_idx-start], &result_batch[batch_idx-start], coord_batch[batch_idx-start].seed_idx, batch_idx);
				}
				else {
					param_batch[batch_idx-start].valid = 0;
				}
			}
			else {
				param_batch[batch_idx-start].valid = 0;
			}
			if (param_batch[batch_idx-start].valid) {
				filled_space += (32 + param_batch[batch_idx-start].string_size);
				*p_readNo = *p_readNo + 1; 
				assert(filled_space <= INPUT_SPACE_SIZE);
				// If the batch is too large to be filled in the space, split them
				/* backup
				if (filled_space + INPUT_THRESHOLD >= INPUT_SPACE_SIZE) {
					if (bwa_verbose >= 4) printf("[debug] 2nd: batch processing on the fly\n");
					// fill input memory space
					fill_input_memory(param_batch, start_read_idx-start, batch_idx-start, input_space, opt);
					// do extension
					extension(input_space, output_space);
					// retreive data from output memory space
					retrieve_output_memory(param_batch, start_read_idx-start, batch_idx-start, output_space, opt);
					// after that, clean the space
					filled_space = 32;
					*p_readNo = 0;
					start_read_idx = batch_idx + 1;
				}
				*/
				if (filled_space + INPUT_THRESHOLD >= INPUT_SPACE_SIZE) {
					if (bwa_verbose >= 4) printf("[debug] 2nd: batch processing on the fly\n");
					if (*p_readNo < MIN_BATCH_SIZE) { // do it by CPU
						int ori_idx;
						for (ori_idx = start_read_idx-start; ori_idx <= batch_idx-start; ori_idx++)
							if (param_batch[ori_idx].valid) original_extension(param_batch[ori_idx].reg, &param_batch[ori_idx], opt);
						// after that, clean the space
						filled_space = 32;
						*p_readNo = 0;
						start_read_idx = batch_idx;
					}
					else {
						// request TBB space for execution
						batch* reservedBatch = reqTBBSpace();
						if (reservedBatch == NULL) { // do it by CPU
							int ori_idx;
							for (ori_idx = start_read_idx-start; ori_idx <= batch_idx-start; ori_idx++)
								if (param_batch[ori_idx].valid) original_extension(param_batch[ori_idx].reg, &param_batch[ori_idx], opt);
							// after that, clean the space
							filled_space = 32;
							*p_readNo = 0;
							start_read_idx = batch_idx;
						}
						else { // do it by FPGA
							// fill input memory space
							int8_t* input_space = (int8_t*)reservedBatch->inputAddr;
							input_space[0] = opt->o_del;
							input_space[1] = opt->e_del;
							input_space[2] = opt->o_ins;
							input_space[3] = opt->e_ins;
							input_space[4] = opt->pen_clip5;
							input_space[5] = opt->pen_clip3;
							input_space[6] = opt->w;
							*((int32_t*)(&input_space[8])) = *p_readNo;
							fill_input_memory(param_batch, start_read_idx-start, batch_idx-start, (int8_t*)reservedBatch->inputAddr, opt);

							// let fpga do the task
							// validate input and send conditional signal
							pthread_mutex_lock(&batchListLock);
							reservedBatch->inputValid = 1;
							pthread_cond_signal(&inputReady);
							//pthread_cond_signal(&reservedBatch->inputReady);
							pthread_mutex_unlock(&batchListLock);

							// retreive data from output memory space
							// wait to get output ready signal
							pthread_mutex_lock(&reservedBatch->batchNodeLock);
							while (!reservedBatch->outputValid) pthread_cond_wait(&reservedBatch->outputReady, &reservedBatch->batchNodeLock);
							pthread_mutex_unlock(&reservedBatch->batchNodeLock);
							retrieve_output_memory(param_batch, (int16_t*)reservedBatch->outputAddr, *p_readNo, start);
							//retrieve_output_memory(param_batch, start_read_idx-start, batch_idx-start, (int8_t*)reservedBatch->outputAddr, opt);
							// after that, clean the space
							releaseBatchSpace(reservedBatch);
							filled_space = 32;
							*p_readNo = 0;
							start_read_idx = batch_idx;
						}
					}
				}
						//original_extension(param_batch[batch_idx-start].reg, &param_batch[batch_idx-start], opt);
			}
		}
		// still has a split to process
		if (bwa_verbose >= 4) printf("[debug] 3rd: batch processing for the tail\n");
		/* back up
		if (filled_space > 32) {
			assert (*p_readNo > 0);
			// fill input memory space
			fill_input_memory(param_batch, start_read_idx-start, batch_idx-1-start, input_space, opt);
			// do extension	
			extension(input_space, output_space);
			// retreive data from output memory space
			retrieve_output_memory(param_batch, start_read_idx-start, batch_idx-1-start, output_space, opt);
			// after that, clean the space
			filled_space = 32;
			*p_readNo = 0;
			start_read_idx = batch_idx;
		}
		*/

		if (filled_space > 32) {
			assert (*p_readNo > 0);
			if (*p_readNo < MIN_BATCH_SIZE) { // do it by CPU
				int ori_idx;
				for (ori_idx = start_read_idx-start; ori_idx <= batch_idx-1-start; ori_idx++)
					if (param_batch[ori_idx].valid) original_extension(param_batch[ori_idx].reg, &param_batch[ori_idx], opt);
				// after that, clean the space
				filled_space = 32;
				*p_readNo = 0;
				start_read_idx = batch_idx;
			}
			else {
				// request TBB space for execution
				batch* reservedBatch = reqTBBSpace();
				if (reservedBatch == NULL) { // do it by CPU
					int ori_idx;
					for (ori_idx = start_read_idx-start; ori_idx <= batch_idx-1-start; ori_idx++)
						if (param_batch[ori_idx].valid) original_extension(param_batch[ori_idx].reg, &param_batch[ori_idx], opt);
					// after that, clean the space
					filled_space = 32;
					*p_readNo = 0;
					start_read_idx = batch_idx;
				}
				else { // do it by FPGA
					// fill input memory space
					int8_t* input_space = (int8_t*)reservedBatch->inputAddr;
					input_space[0] = opt->o_del;
					input_space[1] = opt->e_del;
					input_space[2] = opt->o_ins;
					input_space[3] = opt->e_ins;
					input_space[4] = opt->pen_clip5;
					input_space[5] = opt->pen_clip3;
					input_space[6] = opt->w;
					*((int32_t*)(&input_space[8])) = *p_readNo;
					fill_input_memory(param_batch, start_read_idx-start, batch_idx-1-start, (int8_t*)reservedBatch->inputAddr, opt);

					// let fpga do the task
					// validate input and send conditional signal
					pthread_mutex_lock(&batchListLock);
					reservedBatch->inputValid = 1;
					pthread_cond_signal(&inputReady);
					//pthread_cond_signal(&reservedBatch->inputReady);
					pthread_mutex_unlock(&batchListLock);

					// retreive data from output memory space
					// wait to get output ready signal
					pthread_mutex_lock(&reservedBatch->batchNodeLock);
					while (!reservedBatch->outputValid) pthread_cond_wait(&reservedBatch->outputReady, &reservedBatch->batchNodeLock);
					pthread_mutex_unlock(&reservedBatch->batchNodeLock);
					printf("\n[debug] The output data of the batch is ready, with input address = %llx, and output address = %llx\n", reservedBatch->inputAddr, reservedBatch->outputAddr);

					// FIXME
					retrieve_output_memory(param_batch, (int16_t*)reservedBatch->outputAddr, *p_readNo, start);
					//retrieve_output_memory(param_batch, start_read_idx-start, batch_idx-1-start, (int8_t*)reservedBatch->outputAddr, opt);
					// after that, clean the space
					releaseBatchSpace(reservedBatch);
					filled_space = 32;
					*p_readNo = 0;
					start_read_idx = batch_idx;
				}
			}
		}

		if (bwa_verbose >= 4) printf("[debug] 4th: calculate seed coverage\n");
		// [QA] Calculate Seed Coverage
		for (batch_idx=left_idx; batch_idx<right_idx; batch_idx++) {
			if (param_batch[batch_idx-start].valid)
				mem_chain2aln_seedcov(&chn_batch[batch_idx-start].a[coord_batch[batch_idx-start].chn_idx], param_batch[batch_idx-start].reg);
		}

		if (bwa_verbose >= 4) printf("[debug] 5th: increment coordinates\n");
		// [QA] increment all the coordinates
		for (batch_idx=left_idx; batch_idx<right_idx; batch_idx++) {
			if (chn_batch[batch_idx-start].n > coord_batch[batch_idx-start].chn_idx) {
				if (coord_batch[batch_idx-start].seed_idx > 0) {
					coord_batch[batch_idx-start].seed_idx--;
				}
				else {
					coord_batch[batch_idx-start].chn_idx++;
					coord_batch[batch_idx-start].seed_idx = chn_batch[batch_idx-start].n>coord_batch[batch_idx-start].chn_idx? chn_batch[batch_idx-start].a[coord_batch[batch_idx-start].chn_idx].n-1:  -1;
				}
			}
		}
		if (bwa_verbose >= 4) printf("[debug] 6th: adjust next left and right\n");
		// [QA] Adjust left_idx and right_idx
		while (left_idx < right_idx) {
			if (chn_batch[left_idx-start].n > coord_batch[left_idx-start].chn_idx) {
				break;
			}
			else {
				left_idx++;
			}
		}
		while (left_idx < right_idx) {
			if (chn_batch[right_idx-1-start].n > coord_batch[right_idx-1-start].chn_idx) {
				break;
			}
			else {
				right_idx--;
			}
		}
		if (bwa_verbose >= 4) printf("[debug] complete processing a line of seeds\n");
	}
	if (bwa_verbose >= 4) printf("[debug] all the reads are processed\n");

	for (batch_idx=start; batch_idx<start+batch_size; batch_idx++) {
		for (i = 0; i < chn_batch[batch_idx-start].n; ++i) {
			free(chn_batch[batch_idx-start].a[i].seeds);
		}
		free(chn_batch[batch_idx-start].a);
		regs_batch[batch_idx-start].n = mem_sort_and_dedup(regs_batch[batch_idx-start].n, regs_batch[batch_idx-start].a, opt->mask_level_redun);
		if (opt->flag & MEM_F_NO_EXACT)
			regs_batch[batch_idx-start].n = mem_test_and_remove_exact(opt, regs_batch[batch_idx-start].n, regs_batch[batch_idx-start].a, seqs[batch_idx].l_seq);
		free(param_batch[batch_idx-start].leftQs);
		free(param_batch[batch_idx-start].leftRs);
		free(result_batch[batch_idx-start].srt);
		free(result_batch[batch_idx-start].rseq);
	}
	free(chn_batch);
	free(param_batch);
	free(result_batch);
	free(coord_batch);
	if (bwa_verbose >= 4) printf("[debug] 9th: free and return regs\n");
	return regs_batch;
}

mem_alnreg_v mem_align1_core(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, char *seq)
{
	int i;
	mem_chain_v chn;
	mem_alnreg_v regs;

	for (i = 0; i < l_seq; ++i) // convert to 2-bit encoding if we have not done so
		seq[i] = seq[i] < 4? seq[i] : nst_nt4_table[(int)seq[i]];

	chn = mem_chain(opt, bwt, bns->l_pac, l_seq, (uint8_t*)seq);
	chn.n = mem_chain_flt(opt, chn.n, chn.a);
	if (bwa_verbose >= 4) mem_print_chain(bns, &chn);

	kv_init(regs);
	for (i = 0; i < chn.n; ++i) {
		mem_chain_t *p = &chn.a[i];
		int ret;
		if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
		ret = mem_chain2aln_short(opt, bns->l_pac, pac, l_seq, (uint8_t*)seq, p, &regs);
		if (ret > 0) mem_chain2aln(opt, bns->l_pac, pac, l_seq, (uint8_t*)seq, p, &regs);
		free(chn.a[i].seeds);
	}
	free(chn.a);
	regs.n = mem_sort_and_dedup(regs.n, regs.a, opt->mask_level_redun);
	if (opt->flag & MEM_F_NO_EXACT)
		regs.n = mem_test_and_remove_exact(opt, regs.n, regs.a, l_seq);
	return regs;
}

mem_alnreg_v mem_align1(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq_)
{ // the difference from mem_align1_core() is that this routine: 1) calls mem_mark_primary_se(); 2) does not modify the input sequence
	mem_alnreg_v ar;
	char *seq;
	seq = malloc(l_seq);
	memcpy(seq, seq_, l_seq); // makes a copy of seq_
	ar = mem_align1_core(opt, bwt, bns, pac, l_seq, seq);
	mem_mark_primary_se(opt, ar.n, ar.a, lrand48());
	free(seq);
	return ar;
}

// This routine is only used for the API purpose
mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const char *query_, const mem_alnreg_t *ar)
{
	mem_aln_t a;
	int i, w2, tmp, qb, qe, NM, score, is_rev, last_sc = -(1<<30), l_MD;
	int64_t pos, rb, re;
	uint8_t *query;

	memset(&a, 0, sizeof(mem_aln_t));
	if (ar == 0 || ar->rb < 0 || ar->re < 0) { // generate an unmapped record
		a.rid = -1; a.pos = -1; a.flag |= 0x4;
		return a;
	}
	qb = ar->qb, qe = ar->qe;
	rb = ar->rb, re = ar->re;
	query = malloc(l_query);
	for (i = 0; i < l_query; ++i) // convert to the nt4 encoding
		query[i] = query_[i] < 5? query_[i] : nst_nt4_table[(int)query_[i]];
	a.mapq = ar->secondary < 0? mem_approx_mapq_se(opt, ar) : 0;
	if (ar->secondary >= 0) a.flag |= 0x100; // secondary alignment
	if (bwa_fix_xref2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, opt->w, bns, pac, (uint8_t*)query, &qb, &qe, &rb, &re) < 0) {
		fprintf(stderr, "[E::%s] If you see this message, please let the developer know. Abort. Sorry.\n", __func__);
		exit(1);
	}
	tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
	w2  = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
	w2 = w2 > tmp? w2 : tmp;
	if (bwa_verbose >= 4) printf("* Band width: inferred=%d, cmd_opt=%d, alnreg=%d\n", w2, opt->w, ar->w);
	if (w2 > opt->w) w2 = w2 < ar->w? w2 : ar->w;
//	else w2 = opt->w; // TODO: check if we need this line on long reads. On 1-800bp reads, it does not matter and it should be.
	i = 0; a.cigar = 0;
	do {
		free(a.cigar);
		a.cigar = bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb, (uint8_t*)&query[qb], rb, re, &score, &a.n_cigar, &NM);
		if (bwa_verbose >= 4) printf("* Final alignment: w2=%d, global_sc=%d, local_sc=%d\n", w2, score, ar->truesc);
		if (score == last_sc) break; // it is possible that global alignment and local alignment give different scores
		last_sc = score;
		w2 <<= 1;
	} while (++i < 3 && score < ar->truesc - opt->a);
	l_MD = strlen((char*)(a.cigar + a.n_cigar)) + 1;
	a.NM = NM;
	pos = bns_depos(bns, rb < bns->l_pac? rb : re - 1, &is_rev);
	a.is_rev = is_rev;
	if (a.n_cigar > 0) { // squeeze out leading or trailing deletions
		if ((a.cigar[0]&0xf) == 2) {
			pos += a.cigar[0]>>4;
			--a.n_cigar;
			memmove(a.cigar, a.cigar + 1, a.n_cigar * 4 + l_MD);
		} else if ((a.cigar[a.n_cigar-1]&0xf) == 2) {
			--a.n_cigar;
			memmove(a.cigar + a.n_cigar, a.cigar + a.n_cigar + 1, l_MD); // MD needs to be moved accordingly
		}
	}
	if (qb != 0 || qe != l_query) { // add clipping to CIGAR
		int clip5, clip3;
		clip5 = is_rev? l_query - qe : qb;
		clip3 = is_rev? qb : l_query - qe;
		a.cigar = realloc(a.cigar, 4 * (a.n_cigar + 2) + l_MD);
		if (clip5) {
			memmove(a.cigar+1, a.cigar, a.n_cigar * 4 + l_MD); // make room for 5'-end clipping
			a.cigar[0] = clip5<<4 | 3;
			++a.n_cigar;
		}
		if (clip3) {
			memmove(a.cigar + a.n_cigar + 1, a.cigar + a.n_cigar, l_MD); // make room for 3'-end clipping
			a.cigar[a.n_cigar++] = clip3<<4 | 3;
		}
	}
	a.rid = bns_pos2rid(bns, pos);
	a.pos = pos - bns->anns[a.rid].offset;
	a.score = ar->score; a.sub = ar->sub > ar->csub? ar->sub : ar->csub;
	free(query);
	return a;
}

typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	const mem_pestat_t *pes;
	bseq1_t *seqs;
	mem_alnreg_v *regs;
	int64_t n_processed;
} worker_t;

static void worker1(void *data, int i, int tid)
{
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);
		w->regs[i] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i].l_seq, w->seqs[i].seq);
	} else {
		if (bwa_verbose >= 4) printf("=====> Processing read '%s'/1 <=====\n", w->seqs[i<<1|0].name);
		w->regs[i<<1|0] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i<<1|0].l_seq, w->seqs[i<<1|0].seq);
		if (bwa_verbose >= 4) printf("=====> Processing read '%s'/2 <=====\n", w->seqs[i<<1|1].name);
		w->regs[i<<1|1] = mem_align1_core(w->opt, w->bwt, w->bns, w->pac, w->seqs[i<<1|1].l_seq, w->seqs[i<<1|1].seq);
	}
}

// [QA] The batched-processing version of worker1
static void worker1_batched(void *data, int start, int batch_size, int tid)
{
	worker_t *w = (worker_t*)data;
	mem_alnreg_v* ret = mem_align1_core_batched(w->opt, w->bwt, w->bns, w->pac, w->seqs, start, batch_size);
	int i = -1;
	for (i=start; i<start+batch_size; i++) {
		if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);
		w->regs[i] = ret[i-start];
	}
	free(ret);
}

static void worker2(void *data, int i, int tid)
{
	extern int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2]);
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		if (bwa_verbose >= 4) printf("=====> Finalizing read '%s' <=====\n", w->seqs[i].name);
		mem_mark_primary_se(w->opt, w->regs[i].n, w->regs[i].a, w->n_processed + i);
		mem_reg2sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i], 0, 0);
		free(w->regs[i].a);
	} else {
		if (bwa_verbose >= 4) printf("=====> Finalizing read pair '%s' <=====\n", w->seqs[i<<1|0].name);
		mem_sam_pe(w->opt, w->bns, w->pac, w->pes, (w->n_processed>>1) + i, &w->seqs[i<<1], &w->regs[i<<1]);
		free(w->regs[i<<1|0].a); free(w->regs[i<<1|1].a);
	}
}

void mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	extern void kt_for_batch(int n_threads, void (*func)(void*,int,int,int), void *data, int n, int batch_size); // [QA] new kt_for for batch processing
	worker_t w;
	mem_alnreg_v *regs;
	mem_pestat_t pes[4];
	double ctime, rtime;

	ctime = cputime(); rtime = realtime();
	regs = malloc(n * sizeof(mem_alnreg_v));
	w.opt = opt; w.bwt = bwt; w.bns = bns; w.pac = pac;
	w.seqs = seqs; w.regs = regs; w.n_processed = n_processed;
	w.pes = &pes[0];
	//kt_for(opt->n_threads, worker1, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // find mapping positions
	//[QA] the batched-processing version of worker1(...)
	kt_for_batch(opt->n_threads, worker1_batched, &w, n, opt->batch_size); // find mapping positions
	if (opt->flag&MEM_F_PE) { // infer insert sizes if not provided
		if (pes0) memcpy(pes, pes0, 4 * sizeof(mem_pestat_t)); // if pes0 != NULL, set the insert-size distribution as pes0
		else mem_pestat(opt, bns->l_pac, n, regs, pes); // otherwise, infer the insert size distribution from data
	}
	kt_for(opt->n_threads, worker2, &w, (opt->flag&MEM_F_PE)? n>>1 : n); // generate alignment
	free(regs);
	if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
}
