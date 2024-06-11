#include <memory>
#include <mutex>
#include <cstdlib>   // for getenv
#include <string>    // for std::stoi
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

      // Get the static instance of CublasHandleSingleton class
      static CublasHandleSingleton& GetInstance()
      {
        const char* number_of_openmp_threads = std::getenv("OMP_NUM_THREADS");
        const char* number_of_cuda_devices = std::getenv("CUDA_VISIBLE_DEVICES");
        if (number_of_openmp_threads != nullptr && number_of_cuda_devices != nullptr)
        {
            if (std::stoi(number_of_openmp_threads) > std::stoi(number_of_cuda_devices))
            {
            std::cout << "Mapping multiple OpenMP threads to the same GPU...\n";
            std::lock_guard<std::mutex> lock(GetMutex()); // Lock the mutex so that only one thread will execute the following code at a time;
                                                            // Unlock the mutex when the thread goes out of scope;
                                                            // I think this is needed if multiple threads are mapped to the same GPU.
            }
        }
        else
        {
            throw std::runtime_error("Please set the environment variables OMP_NUM_THREADS and CUDA_VISIBLE_DEVICES.\n");
        }
        static CublasHandleSingleton instance;
        return instance;
      }

      // Method to get the cuBLAS handle
      cublasHandle_t& GetCublasHandle() {
        return handle_;
      }

    private:
      inline static cublasHandle_t handle_;

      static std::mutex& GetMutex()
      {
        static std::mutex mutex;
        return mutex;
      }

      // Private constructor to prevent instantiation
      CublasHandleSingleton()
      {
        // Initialize the cuBLAS handle
        if (!handle_)
        {
          CHECK_CUBLAS_ERROR(cublasCreate(&handle_), "CUBLAS initialization failed...");
        }
      }

      // Private destructor
      ~CublasHandleSingleton()
      {
        // Destroy the cuBLAS handle
        if (handle_)
        {
          CHECK_CUBLAS_ERROR(cublasDestroy(handle_), "CUBLAS finalization failed...");
          this->handle_ = NULL;
        }
      }
  }; // class CublasHandleSingleton
} // namespace micm