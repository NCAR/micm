#include <memory>
#include <mutex>
#include <map>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <micm/util/cuda_util.cuh>

namespace micm
{
  class CublasHandleSingleton {
    public:
      // Delete copy constructor and assignment operator to prevent copies
      CublasHandleSingleton(const CublasHandleSingleton&) = delete;
      CublasHandleSingleton& operator=(const CublasHandleSingleton&) = delete;

      // Get the static instance of CublasHandleSingleton class
      static CublasHandleSingleton& GetInstance()
      {
        int device_id;
        CHECK_CUDA_ERROR(cudaGetDevice(&device_id), "Failed to get device ID...");
        if (auto search = cublas_handle_map_.find(device_id); search == cublas_handle_map_.end())
        {
          std::lock_guard<std::mutex> lock(GetMutex()); // no cublas handle if found; lock the mutex and generate a new cublas handle below
        }
        static CublasHandleSingleton instance;          // create the cublas handle inside
        if (auto search = cublas_handle_map_.find(device_id); search == cublas_handle_map_.end())
        {
          cublas_handle_map_[device_id] = handle_;      // save the cublas handle to the map
        }
        return instance;
      }

      // Method to get the cuBLAS handle
      cublasHandle_t& GetCublasHandle()
      {
        return handle_;
      }

    private:
      inline static cublasHandle_t handle_;
      inline static std::map<int, cublasHandle_t> cublas_handle_map_;

      static std::mutex& GetMutex()
      {
        static std::mutex mutex;
        return mutex;
      }

      // Private constructor to prevent instantiation
      CublasHandleSingleton()
      {
        // Initialize the cublas handle map
        if (cublas_handle_map_.empty())
        {
          cublas_handle_map_ = {};
        }
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
        // Destroy the cublas handle map
        if (!cublas_handle_map_.empty())
        {
          cublas_handle_map_.clear();
        }
      }
  }; // class CublasHandleSingleton
} // namespace micm