#ifdef HAVE_CONFIG_H
# include <config.h>
#endif // HAVE_CONFIG_H
#include <stdio.h>

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

#ifdef HAVE_STRING_H
# include <string.h>
#endif // HAVE_STRING_H

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif // HAVE_UNISTD_H

#include <pthread.h>
#include <string.h>
#include <iostream>
#include <iomanip>
// #include <fstream>
#include <algorithm>

#include <aalsdk/ccilib/CCILib.h>
#include <aalsdk/aalclp/aalclp.h>

USING_NAMESPACE(std)
USING_NAMESPACE(AAL)
USING_NAMESPACE(CCILib)

#ifndef CL
# define CL(x)                     ((x) * 64)
#endif // CL
#ifndef LOG2_CL
# define LOG2_CL                   6
#endif // LOG2_CL
#ifndef MB
# define MB(x)                     ((x) * 1024 * 1024)
#endif // MB

#define CACHELINE_ALIGNED_ADDR(p)  ((p) >> LOG2_CL)
#define LPBK1_BUFFER_SIZE          CL(1)

#define BWA_NUM_BATCHES            4
#define BWA_INPUT_BUFFER_SIZE      CL(4096)     // the size of TBB
#define BWA_OUTPUT_BUFFER_SIZE     CL(256)      // the size of RBB

#define LPBK1_DSM_SIZE             MB(4)

#define CSR_CIPUCTL                0x280

#define CSR_AFU_DSM_BASEL          0x1a00
#define CSR_AFU_DSM_BASEH          0x1a04
#define CSR_SRC_ADDR               0x1a20
#define CSR_DST_ADDR               0x1a24
#define CSR_NUM_LINES              0x1a28
#define CSR_CTL                    0x1a2c
#define CSR_CFG                    0x1a34

#define CSR_OFFSET(x)              ((x) / sizeof(bt32bitCSR))

#define DSM_STATUS_TEST_COMPLETE   0x40
#define DSM_STATUS_TEST_ERROR      0x44
#define DSM_STATUS_MODE_ERROR_0    0x60

#define DSM_STATUS_ERROR_REGS      8

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

#define NUM_THREADS 12

pthread_mutex_t batchLock = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;  // double lock possibility
pthread_cond_t  gotBatch  = PTHREAD_COND_INITIALIZER;

int numBatches = 0;

typedef struct thdata
{
   int thread_id;
   ICCIDevice* device;
} thdata;

struct batch {
    btVirtAddr batchAddr;
    struct batch* next;
};

struct batch* batchList = NULL;
struct batch* lastBatch = NULL;

void 
addBatch(btVirtAddr batchAddr, 
        pthread_mutex_t* p_mutex, 
        pthread_cond_t* p_cond_var)
{
    int rc;
    struct batch* bat;

    bat = new struct batch;
    bat->batchAddr = batchAddr;
    bat->next = NULL;

    // Lock the mutex, to assure exclusive access to the list
    rc = pthread_mutex_lock(p_mutex);

    if (numBatches == 0) {  // list is empty
        batchList = bat;
        lastBatch = bat;
    } else {
        lastBatch->next = bat;
        lastBatch = bat;
    }

    numBatches++;

    // Unlock mutex
    rc = pthread_mutex_unlock(p_mutex);
    // Signal the condition variable - there's a new batch available
    rc = pthread_cond_signal(p_cond_var);
}

struct batch* 
getBatch(pthread_mutex_t* p_mutex)
{
    int rc;
    struct batch* bat;

    // Lock the mutex, to assure exclusive access to the list
    rc = pthread_mutex_lock(p_mutex);

    if (numBatches > 0) {
        bat = batchList;
        batchList = bat->next;
        if (batchList == NULL) {    // This was the last batch on the list
            lastBatch = NULL;
        }
        numBatches--;
    } else {    // Batch list is empty
        bat = NULL;
    }

    // Unlock mutex
    rc = pthread_mutex_unlock(p_mutex);

    return bat;
}

void
worker(struct batch* bat, int thread_id)
{
    if (bat) {
        cout << "Thread " << thread_id << " acquired batch in "
            << bat->batchAddr << endl;
    }
}

void*
loop_worker(void* data)
{
    int rc;
    struct batch* bat;
    int thread_id = *((int*)data);

    rc = pthread_mutex_lock(&batchLock);

    while (true) {
        if (numBatches > 0) {
            bat = getBatch(&batchLock);
            if (bat) {
                worker(bat, thread_id);
                addBatch(bat->batchAddr, &batchLock, &gotBatch);
                free(bat);
            }
        } else {
            rc = pthread_cond_wait(&gotBatch, &batchLock);
        }
    }
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

    ICCIWorkspace *pDSMWorkspace    = pCCIDevice->AllocateWorkspace(LPBK1_DSM_SIZE);
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

    // Assert Device Reset
    pCCIDevice->SetCSR(CSR_CTL, 0);

    // Clear the DSM
    memset((void *)pDSMUsrVirt, 0, pDSMWorkspace->GetSizeInBytes());

    // De-assert Device Reset
    pCCIDevice->SetCSR(CSR_CTL, 1);

    // Set input workspace address
    pCCIDevice->SetCSR(CSR_SRC_ADDR, CACHELINE_ALIGNED_ADDR(pInputWorkspace->GetPhysicalAddress()));

    // Set output workspace address
    pCCIDevice->SetCSR(CSR_DST_ADDR,  CACHELINE_ALIGNED_ADDR(pOutputWorkspace->GetPhysicalAddress()));

    // Set the test mode
    pCCIDevice->SetCSR(CSR_CFG,       0);

    volatile bt32bitCSR *StatusAddr = (volatile bt32bitCSR *)
                                    (pDSMUsrVirt  + DSM_STATUS_TEST_COMPLETE);

    // Start the test
    pCCIDevice->SetCSR(CSR_CTL,      0x3);

    int       thread_id[NUM_THREADS];
    pthread_t pThreads[NUM_THREADS];
    for (i = 0; i < NUM_THREADS; ++i) {
        thread_id[i] = i;
        pthread_create(&pThreads[i], NULL, loop_worker, (void*)&thread_id[i]);
    }
    sleep(3);

    addBatch(pInputUsrVirt, &batchLock, &gotBatch);
    addBatch(pInputUsrVirt + BWA_INPUT_BUFFER_SIZE / sizeof(btUnsigned32bitInt), 
      &batchLock, &gotBatch);
    addBatch(pInputUsrVirt + BWA_INPUT_BUFFER_SIZE * 2 / sizeof(btUnsigned32bitInt), 
      &batchLock, &gotBatch);
    addBatch(pInputUsrVirt + BWA_INPUT_BUFFER_SIZE * 3 / sizeof(btUnsigned32bitInt), 
      &batchLock, &gotBatch);

    // Wait for test completion
    // while( BWA_NUM_BATCHES != *StatusAddr ) {
    //     usleep(100);
    // }

    sleep(1);

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
