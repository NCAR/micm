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
      static CublasHandleSingleton& GetInstance(int device_id)
      {
        if (auto search = cublas_handle_map_.find(device_id); search == cublas_handle_map_.end())
        {
          std::lock_guard<std::mutex> lock(GetMutex()); // no cublas handle if found; lock the mutex and generate a new cublas handle below
        }
        static CublasHandleSingleton instance;
        if (auto search = cublas_handle_map_.find(device_id); search == cublas_handle_map_.end())
        {
          cublasHandle_t handle;
          CHECK_CUBLAS_ERROR(cublasCreate(&handle), "CUBLAS initialization failed..."); // create the cublas handle
          cublas_handle_map_[device_id] = handle;                                       // save the cublas handle
        }
        return instance;
      }

      // Method to get the cuBLAS handle
      cublasHandle_t& GetCublasHandle(int device_id)
      {
        return cublas_handle_map_[device_id];
      }

    private:
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
      }

      // Private destructor
      ~CublasHandleSingleton()
      {
        // Destroy the cublas handle map
        if (!cublas_handle_map_.empty())
        {
          for (const auto& pair : cublas_handle_map_)
          {
            CHECK_CUBLAS_ERROR(cublasDestroy(pair.second), "CUBLAS finalization failed..."); // destroy the cublas handle
          }
          cublas_handle_map_.clear();
        }
      }
  }; // class CublasHandleSingleton
} // namespace micm