#include <iostream> 
struct pair{
    size_t first; 
    size_t second; 
}
namespace micm{
    namespace cuda{
        __global__ void Decompose_kernel(){

        }
    
        void decompose_kernelSetup(
            const double* A, 
            double* L, 
            double* U,
            size_t* niLu, //pair 
            bool* do_aik,
            size_t* aik, 
            size_t* uik_nkj, //pair 
            size_t* lij_ujk, //pair 
            bool* do_aki, 
            size_t* aki, 
            size_t* lki_nkj, //pair 
            size_t* lkj_uji, //pair
            size_t* uii){
            
            vector<


            }
        


    }
}