#include <memory>
#include <mutex>
#include <cublas_v2.h>
#include <micm/util/cuda_util.cuh>

namespace micm
{
    class CublasHandleSingleton {
public:
    // Delete copy constructor and assignment operator to prevent copies
    CublasHandleSingleton(const CublasHandleSingleton&) = delete;
    CublasHandleSingleton& operator=(const CublasHandleSingleton&) = delete;

    // Static method to get the singleton instance
    static CublasHandleSingleton& GetInstance() {
        static CublasHandleSingleton instance;
        return instance;
    }

    // Method to get the cublas handle
    cublasHandle_t & GetCublasHandle() {
        return handle;
    }

private:
    inline static cublasHandle_t handle;

    // Private constructor
    CublasHandleSingleton() {
        // Initialize the cuBLAS handle
        std::cout << "JS: ready to call constructor..." << std::endl;
        cublasStatus_t status = cublasCreate(&handle);
        if (status != CUBLAS_STATUS_SUCCESS) {
            std::cerr << "CUBLAS initialization failed\n";
            throw std::runtime_error("CUBLAS initialization failed");
        }
    }

    // Private destructor
    ~CublasHandleSingleton() {
        // Destroy the cuBLAS handle
        std::cout << "JS: ready to call destructor..." << std::endl;
        cublasDestroy(handle);
    }
};


//     class CublasHandleSingleton {
//     public:
//         // Get the static instance of CublasHandleSingleton class
//         static CublasHandleSingleton& GetInstance() {
//             std::lock_guard<std::mutex> lock(GetMutex()); // Lock the mutex
//             std::cout << "JS: ready to call constructor..." << std::endl;
//             static CublasHandleSingleton instance;
//             return instance;
//         }

//         // Get the cublas handle as a shared pointer with custom deleter
//         cublasHandle_t GetCublasHandle() {
//             std::cout << "JS: GetCublasHandle is called" << std::endl;
//             return handle_;
//         }

//     private:
//         cublasHandle_t handle_;
//         static std::mutex& GetMutex() {
//             static std::mutex mutex;
//             return mutex;
//         }

//         // Private constructor to prevent instantiation
//         CublasHandleSingleton()
//         {
// //            if (!handle_)
// //            {
//                 CHECK_CUBLAS_ERROR(cublasCreate(&handle_), "CUBLAS initialization failed...");
// //            }
//         }
        
//         // Delete copy constructor and assignment operator
//         CublasHandleSingleton(const CublasHandleSingleton&) = delete;
//         CublasHandleSingleton& operator=(const CublasHandleSingleton&) = delete;

//         // Destructor
//         ~CublasHandleSingleton()
//         {
//             if (handle_)
//             {
//                 CHECK_CUBLAS_ERROR(cublasDestroy(handle_), "CUBLAS finalization failed...");
//                 this->handle_ = NULL;
//             }
//         }
//     };
}