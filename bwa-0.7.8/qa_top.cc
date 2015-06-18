#ifdef HAVE_CONFIG_H
# include <config.h>
#endif // HAVE_CONFIG_H

# include <config.h>

#include <stdio.h>
#include <assert.h>

// Real one, comment for testing
#ifdef STDC_HEADERS
# include <stdlib.h>
# include <stddef.h>
#else
# ifdef HAVE_STDLIB_H
#    include <stdlib.h>
# else
#    error Required system header stdlib.h not found.
# endif // HAVE_STDLIB_H
#endif // STDC_HEADERS

//just for testing, should be deleted
// # include <stdlib.h>
// # include <stddef.h>
// # include <config.h>


#ifdef HAVE_STRING_H
# include <string.h>
#endif // HAVE_STRING_H

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif // HAVE_UNISTD_H

#include <pthread.h>

#include <iostream>
#include <iomanip>
// #include <fstream>
#include <algorithm>

#include "qa_top.h"

#include <aalsdk/ccilib/CCILib.h>
#include <aalsdk/aalclp/aalclp.h>

USING_NAMESPACE(std)
USING_NAMESPACE(AAL)
USING_NAMESPACE(CCILib)

extern "C" {
#include "batch.h"
	void extension(int8_t* input_space, int8_t* output_space);
	int bwa_main(int argc, char* argv[]);
	batch* reqTBBSpace(); 
	void releaseBatchSpace(batch*);
}

///////////////////////////////////////////////
BEGIN_C_DECLS

struct CCIDemoCmdLine
{
    btUIntPtr               flags;
#define CCIDEMO_CMD_FLAG_HELP       0x00000001
#define CCIDEMO_CMD_FLAG_VERSION    0x00000002

    CCIDeviceImplementation target;
    int                     log;
    int                     trace;
};

struct CCIDemoCmdLine gCCIDemoCmdLine = 
{
    0,
    CCI_NULL,
    0,
    0
};

int ccidemo_on_nix_long_option_only(AALCLP_USER_DEFINED , const char * );
int ccidemo_on_nix_long_option(AALCLP_USER_DEFINED , const char * , const char * );

aalclp_option_only ccidemo_nix_long_option_only = { ccidemo_on_nix_long_option_only, };
aalclp_option      ccidemo_nix_long_option      = { ccidemo_on_nix_long_option,      };

void help_msg_callback(FILE * , struct _aalclp_gcs_compliance_data * );
void showhelp(FILE * , struct _aalclp_gcs_compliance_data * );

AALCLP_DECLARE_GCS_COMPLIANT(stdout,
                             "CCIDemo",
                             CCILIB_VERSION,
                             "",
                             help_msg_callback,
                             &gCCIDemoCmdLine)

int parsecmds(struct CCIDemoCmdLine * , int , char *[] );
int verifycmds(struct CCIDemoCmdLine * );

END_C_DECLS
///////////////////////////////////////////////


// global variables
// ? Do we need recursive lock
//pthread_mutex_t batchLock = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;  // double lock possibility

pthread_mutex_t batchListLock = PTHREAD_MUTEX_INITIALIZER;  // standard mutex initializer
pthread_mutex_t freeListLock = PTHREAD_MUTEX_INITIALIZER;  // standard mutex initializer
pthread_cond_t inputReady = PTHREAD_COND_INITIALIZER;

// initialize the linked list for occupied batches and free batches
batch* batchListDummyHead = NULL;
batch* batchListTail = NULL;
batch* freeListDummyHead = NULL;
batch* freeListTail = NULL;
batch  bat[BWA_NUM_BATCHES];

int request_pearray[BWA_NUM_BATCHES];

volatile bt32bitCSR *StatusAddr;

// end global

batch* 
reqTBBSpace()
{
    int rc;
    batch* p_newBatch = NULL;
    batch* p_nextBatch = NULL;

    // Lock the free list, find if there is at least one free node available
    rc = pthread_mutex_lock(&freeListLock);

    if (freeListDummyHead->next == NULL) {
        // no free node, then return null pointer to let CPU do the job
        // 1st: unlock the free list
        // 2nd: return NULL
        rc = pthread_mutex_unlock(&freeListLock);
	return NULL;
    }
    else {
        // there is at least one free node, take it for the requesting thread

        // 1st: remove one node from the free list
        p_newBatch = freeListDummyHead->next;
	// cond 1: there is only one node, meaning the tail should be modified
	if (p_newBatch == freeListTail) {
		freeListDummyHead->next = NULL;
		freeListTail = freeListDummyHead;
	}
	//cond 2: there is not only one node, meaning the tail should not be modified
	else {
		p_nextBatch = p_newBatch->next;
        	freeListDummyHead->next = p_nextBatch;
		p_nextBatch->prev = freeListDummyHead;
	}
	// unlock the free list, this finishes the first step
	rc = pthread_mutex_unlock(&freeListLock);

        // 2nd: add the free node to the tail of the batch list
        // this step needs to lock the batch list
        rc = pthread_mutex_lock(&batchListLock);
	batchListTail->next = p_newBatch;
	p_newBatch->prev = batchListTail;
	p_newBatch->next = NULL;
	batchListTail = p_newBatch;
	rc = pthread_mutex_unlock(&batchListLock);
    }

    return p_newBatch;
}

void releaseBatchSpace(batch* p_batch) {
    int rc;

    // for release, add the batch node back to the free list
    pthread_mutex_t* p_mutex = &freeListLock;

    p_batch->inputValid = 0;
    p_batch->outputValid = 0;

    // Lock the mutex, to assure exclusive access to the list
    rc = pthread_mutex_lock(p_mutex);

    // add to the tail of the free list
    // it has already been removed from the batch list during execution
    freeListTail->next = p_batch;
    p_batch->prev = freeListTail;
    p_batch->next = NULL;
    freeListTail = p_batch;

    rc = pthread_mutex_unlock(p_mutex);
}

void show_req_vector()
{
    cout << "    Request Vector = ";
    for (int k = 0; k < BWA_NUM_BATCHES; k++) 
        cout << request_pearray[k];
    cout << endl;
}

void show_status(bt32bitCSR status)
{
    cout << "Current Status = ";
    for (int k = 0; k < BWA_NUM_BATCHES; k++)
        cout << ((status >> k) & 1);
    cout << endl;
}

void*
execution_fpga(void* data)
{
    // char* exe_message = (char*) data;
    //printf ("Start running the execution_fpga thread with message: %s\n", exe_message);

    cout << "Start running the execution thread." << endl;

    ICCIDevice *pCCIDevice = (ICCIDevice*) data;

    int rc;
    batch* cur_batch = NULL;
    batch* next_batch = NULL;

    for (int k = 0; k < BWA_NUM_BATCHES; k++)
        request_pearray[k] = 0;

    while (1) {
        // 1st: wait for a valid input
        rc = pthread_mutex_lock(&batchListLock);
        while (batchListDummyHead == NULL || batchListDummyHead->next == NULL || batchListDummyHead->next->inputValid == 0) {
            rc = pthread_cond_wait(&inputReady, &batchListLock);
        }

        // 2nd: remove the node from the batch list
        cur_batch = batchListDummyHead->next;
        assert (cur_batch->inputValid && !cur_batch->outputValid);
        // cond 1: the current batch is the tail
        if (cur_batch == batchListTail) {
            //then, remove the only node from the batch list
            batchListDummyHead->next = NULL;
            batchListTail = batchListDummyHead;
        }
        // cond 2: the current batch is not the tail
        else {
            //then, remove the current node
            next_batch = cur_batch->next;
            batchListDummyHead->next = next_batch;
            next_batch->prev = batchListDummyHead;
        }
        rc = pthread_mutex_unlock(&batchListLock);

        cout << "[" << cur_batch->idx << "] ";
        cout << "Batch (" << *(btUnsigned32bitInt *)cur_batch->inputAddr << ") dispatched to FPGA" << endl;

        request_pearray[cur_batch->idx] = 1;

        int request_vector = 0;
        for (int k = 0; k < BWA_NUM_BATCHES; k++)
            request_vector |= ( request_pearray[k] << k );
        show_req_vector();

        // pCCIDevice->SetCSR(CSR_NUM_LINES, 1);
        pCCIDevice->SetCSR(CSR_REQ_PEARRAY, request_vector);
    }
    
    return NULL;
}

void*
collection_fpga(void* data)
{
    cout << "Start running the collection thread." << endl;

    ICCIDevice *pCCIDevice = (ICCIDevice*) data;

    bt32bitCSR pre_status = 0;
    int pre_pearray_busy[BWA_NUM_BATCHES];
    int cur_pearray_busy[BWA_NUM_BATCHES];
    bool request_updtd = true;
    int rc;

    for (int k = 0; k < BWA_NUM_BATCHES; k++) {
        pre_pearray_busy[k] = 0;
        cur_pearray_busy[k] = 0;
    }

    while (true) {
        while ( pre_status == *StatusAddr ) {
            usleep(100);
        }

        pre_status = *StatusAddr;

        cout << endl << "FPGA STATUS CHANGED ";
        show_status(pre_status);

        for (int k = 0; k < BWA_NUM_BATCHES; k++) {
            cur_pearray_busy[k] = (pre_status >> k) & 1;
            // cout << cur_pearray_busy[k];

            if (pre_pearray_busy[k] == 0 && cur_pearray_busy[k] == 1) {
                request_pearray[k] = 0;     // no race condition with the exe thread
                request_updtd = false;

                cout << "[" << k << "] ";
                cout << "PEarray starts to work. Updating request vector." << endl;
            } else if (pre_pearray_busy[k] == 1 && cur_pearray_busy[k] == 0) {
                rc = pthread_mutex_lock(&bat[k].batchNodeLock);
                bat[k].outputValid = 1;
                pthread_cond_signal(&bat[k].outputReady);
                rc = pthread_mutex_unlock(&bat[k].batchNodeLock);

                cout << "[" << k << "] ";
                cout << "PEarray finished processing. Releasing current batch." << endl;
            }

            pre_pearray_busy[k] = cur_pearray_busy[k];
        }

        if (!request_updtd) {
            request_updtd = true;
            int request_vector = 0;
            for (int k = 0; k < BWA_NUM_BATCHES; k++)
                request_vector |= ( request_pearray[k] << k );
            show_req_vector();

            // update the request vector
            pCCIDevice->SetCSR(CSR_REQ_PEARRAY, request_vector);
        }
    }

    return NULL;
}

int main(int argc, char *argv[])
{

    if (argc < 2) {
        showhelp(stdout, &_aalclp_gcs_data);
        return 1;
    } else if (parsecmds(&gCCIDemoCmdLine, argc, argv)) {
        cerr << "Error scanning command line." << endl;
        return 2;
    } else if (flag_is_set(gCCIDemoCmdLine.flags, CCIDEMO_CMD_FLAG_HELP | CCIDEMO_CMD_FLAG_VERSION)) {
        return 0;
    } else if (verifycmds(&gCCIDemoCmdLine)) {
        return 3;
    }

    const CCIDeviceImplementation CCIDevImpl = gCCIDemoCmdLine.target;

    ICCIDeviceFactory *pCCIDevFactory = GetCCIDeviceFactory(CCIDevImpl);

    ICCIDevice *pCCIDevice = pCCIDevFactory->CreateCCIDevice();

#if (1 == ENABLE_DEBUG)
    pCCIDevice->GetSynchronizer()->SetLogLevel(gCCIDemoCmdLine.log);
    pCCIDevice->GetSynchronizer()->SetTraceLevel(gCCIDemoCmdLine.trace);
#endif // ENABLE_DEBUG

    ICCIWorkspace *pDSMWorkspace    = pCCIDevice->AllocateWorkspace(BWA_DSM_SIZE);
    ICCIWorkspace *pInputWorkspace  = pCCIDevice->AllocateWorkspace(BWA_INPUT_BUFFER_SIZE * BWA_NUM_BATCHES);
    ICCIWorkspace *pOutputWorkspace = pCCIDevice->AllocateWorkspace(BWA_OUTPUT_BUFFER_SIZE * BWA_NUM_BATCHES);

    volatile btVirtAddr pInputUsrVirt  = pInputWorkspace->GetUserVirtualAddress(); 
    volatile btVirtAddr pOutputUsrVirt = pOutputWorkspace->GetUserVirtualAddress();
    volatile btVirtAddr pDSMUsrVirt    = pDSMWorkspace->GetUserVirtualAddress();

    memset((void *)pOutputUsrVirt, 0, pOutputWorkspace->GetSizeInBytes());
    memset((void *)pDSMUsrVirt, 0, pDSMWorkspace->GetSizeInBytes());

    bt32bitCSR i;
    bt32bitCSR csr;

    // Assert CAFU Reset
    csr = 0;
    pCCIDevice->GetCSR(CSR_CIPUCTL, &csr);
    csr |= 0x01000000;
    pCCIDevice->SetCSR(CSR_CIPUCTL, csr);

    // De-assert CAFU Reset
    csr = 0;
    pCCIDevice->GetCSR(CSR_CIPUCTL, &csr);
    csr &= ~0x01000000;
    pCCIDevice->SetCSR(CSR_CIPUCTL, csr);

    // Set DSM base, high then low
    pCCIDevice->SetCSR(CSR_AFU_DSM_BASEH, pDSMWorkspace->GetPhysicalAddress() >> 32);
    pCCIDevice->SetCSR(CSR_AFU_DSM_BASEL, pDSMWorkspace->GetPhysicalAddress() & 0xffffffff);

    // Poll for AFU ID
    do
    {
        csr = *(volatile btUnsigned32bitInt *)pDSMUsrVirt;
    } while( 0 == csr );

    // Print the AFU ID
    cout << "AFU ID=";
    for ( i = 0 ; i < 4 ; ++i ) {
        cout << std::setw(8) << std::hex << std::setfill('0')
            << *(btUnsigned32bitInt *)(pDSMUsrVirt + (3 - i) * sizeof(btUnsigned32bitInt));
    }
    cout << endl;

    // Assert Device Reset
    pCCIDevice->SetCSR(CSR_CTL, 0);

    // Clear the DSM
    memset((void *)pDSMUsrVirt, 0, pDSMWorkspace->GetSizeInBytes());

    // De-assert Device Reset
    pCCIDevice->SetCSR(CSR_CTL, 1);

    // Set input workspace address
    pCCIDevice->SetCSR(CSR_SRC_ADDR, CACHELINE_ALIGNED_ADDR(pInputWorkspace->GetPhysicalAddress()));

    // Set output workspace address
    pCCIDevice->SetCSR(CSR_DST_ADDR, CACHELINE_ALIGNED_ADDR(pOutputWorkspace->GetPhysicalAddress()));

    // Set the test mode
    pCCIDevice->SetCSR(CSR_CFG, 0);

    StatusAddr = (volatile bt32bitCSR *)(pDSMUsrVirt  + DSM_STATUS_PEARRAY);

    pCCIDevice->SetCSR(CSR_REQ_PEARRAY, 0);

    // Start the test
    pCCIDevice->SetCSR(CSR_CTL, 0x3);


    // [QA] Initialize the batch list and the free list
    batchListDummyHead = new batch;
    batchListDummyHead->next = NULL;
    batchListDummyHead->prev = NULL;
    batchListTail = batchListDummyHead;

    freeListDummyHead = new batch;
    freeListDummyHead->next = NULL;
    freeListDummyHead->prev = NULL;
    freeListTail = freeListDummyHead;

    for (int k = 0; k < BWA_NUM_BATCHES; k++) {
        freeListTail->next = &bat[k];
        freeListTail->next->prev = freeListTail;
        freeListTail = freeListTail->next;
        freeListTail->inputValid = 0;
        freeListTail->outputValid = 0;
        freeListTail->idx = k;
        
        freeListTail->inputAddr = pInputUsrVirt + BWA_INPUT_BUFFER_SIZE * k;  // byte addressing
        freeListTail->outputAddr = pOutputUsrVirt + BWA_OUTPUT_BUFFER_SIZE * k;

        freeListTail->next = NULL;
        //freeListTail->inputReady = PTHREAD_COND_INITIALIZER;
        freeListTail->outputReady = PTHREAD_COND_INITIALIZER;
        freeListTail->batchNodeLock = PTHREAD_MUTEX_INITIALIZER;
    }

    cout << "done" << endl;

    pthread_t exe_thread;
    pthread_t col_thread;
    pthread_create(&exe_thread, NULL, execution_fpga,  (void*) pCCIDevice);
    pthread_create(&col_thread, NULL, collection_fpga, (void*) pCCIDevice);

    /* skip the first command */
    bwa_main(argc - 1, argv + 1);

    batch* p = NULL;
    while (batchListDummyHead) {
        p = batchListDummyHead;
        batchListDummyHead = batchListDummyHead->next;
        delete p;
    }
    while (freeListDummyHead) {
        p = freeListDummyHead;
        freeListDummyHead = freeListDummyHead->next;
        delete p;
    }

    // Stop the device
    pCCIDevice->SetCSR(CSR_CTL,      0x7);

    // Release the Workspaces
    pCCIDevice->FreeWorkspace(pInputWorkspace);
    pCCIDevice->FreeWorkspace(pOutputWorkspace);
    pCCIDevice->FreeWorkspace(pDSMWorkspace);

    // Release the CCI Device instance.
    pCCIDevFactory->DestroyCCIDevice(pCCIDevice);

    // Release the CCI Device Factory instance.
    delete pCCIDevFactory;

    return 0;
}


BEGIN_C_DECLS

int ccidemo_on_nix_long_option_only(AALCLP_USER_DEFINED user, const char *option)
{
   struct CCIDemoCmdLine *cl = (struct CCIDemoCmdLine *)user;

   if ( 0 == strcmp("--help", option) ) {
      flag_setf(cl->flags, CCIDEMO_CMD_FLAG_HELP);
   } else if ( 0 == strcmp("--version", option) ) {
      flag_setf(cl->flags, CCIDEMO_CMD_FLAG_VERSION);
   }

   return 0;
}

int ccidemo_on_nix_long_option(AALCLP_USER_DEFINED user, const char *option, const char *value)
{
   struct CCIDemoCmdLine *cl = (struct CCIDemoCmdLine *)user;

   if ( 0 == strcmp("--target", option) ) {
      if ( 0 == strcasecmp("aal", value) ) {
#if (1 == CCILIB_ENABLE_AAL)
         cl->target = CCI_AAL;
#else
         cout << "The version of CCILib was built without support for --target=AAL" << endl;
         return 1;
#endif // CCILIB_ENABLE_AAL
      } else if ( 0 == strcasecmp("ase", value) ) {
#if (1 == CCILIB_ENABLE_ASE)
         cl->target = CCI_ASE;
#else
         cout << "The version of CCILib was built without support for --target=ASE" << endl;
         return 2;
#endif // CCILIB_ENABLE_ASE
      } else if ( 0 == strcasecmp("direct", value) ) {
#if (1 == CCILIB_ENABLE_DIRECT)
         cl->target = CCI_DIRECT;
#else
         cout << "The version of CCILib was built without support for --target=Direct" << endl;
         return 3;
#endif // CCILIB_ENABLE_DIRECT
      } else {
         cout << "Invalid value for --target : " << value << endl;
         return 4;
      }
   } else if ( 0 == strcmp("--log", option) ) {
      char *endptr = NULL;

      cl->log = (int)strtol(value, &endptr, 0);
      if ( endptr != value + strlen(value) ) {
         cl->log = 0;
      } else if ( cl->log < 0) {
         cl->log = 0;
      } else if ( cl->log > 7) {
         cl->log = 7;
      }
   } else if ( 0 == strcmp("--trace", option) ) {
      char *endptr = NULL;

      cl->trace = (int)strtol(value, &endptr, 0);
      if ( endptr != value + strlen(value) ) {
         cl->trace = 0;
      } else if ( cl->trace < 0) {
         cl->trace = 0;
      } else if ( cl->trace > 7) {
         cl->trace = 7;
      }
   }

   return 0;
}

void help_msg_callback(FILE *fp, struct _aalclp_gcs_compliance_data *gcs)
{
   fprintf(fp, "Usage:\n");
   fprintf(fp, "   CCIDemo [--target=<TARGET>]\n");
   fprintf(fp, "\n");
   fprintf(fp, "      <TARGET> = one of { ");
#if (1 == CCILIB_ENABLE_AAL)
   fprintf(fp, "AAL ");
#endif // CCILIB_ENABLE_AAL
#if (1 == CCILIB_ENABLE_ASE)
   fprintf(fp, "ASE ");
#endif // CCILIB_ENABLE_ASE
#if (1 == CCILIB_ENABLE_DIRECT)
   fprintf(fp, "Direct ");
#endif // CCILIB_ENABLE_DIRECT
   fprintf(fp, "}\n");
   fprintf(fp, "\n");
}

void showhelp(FILE *fp, struct _aalclp_gcs_compliance_data *gcs)
{
   help_msg_callback(fp, gcs);
}

int parsecmds(struct CCIDemoCmdLine *cl, int argc, char *argv[])
{
   int    res;
   int    clean;
   aalclp clp;

   res = aalclp_init(&clp);
   if ( 0 != res ) {
      cerr << "aalclp_init() failed : " << res << ' ' << strerror(res) << endl;
      return res;
   }

   ccidemo_nix_long_option_only.user = cl;
   aalclp_add_nix_long_option_only(&clp, &ccidemo_nix_long_option_only);

   ccidemo_nix_long_option.user = cl;
   aalclp_add_nix_long_option(&clp, &ccidemo_nix_long_option);

   res = aalclp_add_gcs_compliance(&clp);
   if ( 0 != res ) {
      cerr << "aalclp_add_gcs_compliance() failed : " << res << ' ' << strerror(res) << endl;
      goto CLEANUP;
   }

   res = aalclp_scan_argv(&clp, argc, argv);
   if ( 0 != res ) {
      cerr << "aalclp_scan_argv() failed : " << res << ' ' << strerror(res) << endl;
   }

CLEANUP:
   clean = aalclp_destroy(&clp);
   if ( 0 != clean ) {
      cerr << "aalclp_destroy() failed : " << clean << ' ' << strerror(clean) << endl;
   }

   return res;
}

int verifycmds(struct CCIDemoCmdLine *cl)
{
   if ( CCI_NULL == cl->target ) {
      cout << "No valid --target specified." << endl;
      return 1;
   }

   return 0;
}

END_C_DECLS

