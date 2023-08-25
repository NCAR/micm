#pragma once
#include <vector> 
namespace micm{
    class CUDAMatrixParam{
        public: 
        const double* rate_constants_; 
        const double* state_variables_; 
        double* forcing_; 
        double* jacobian_; 
        size_t n_grids_; 
        size_t n_reactions_; 
        size_t n_species_; 
        size_t jacobian_size_; 
        
        CUDAMatrixParam(){
        
        }; 
        inline void setGrids(size_t n_grids){
            n_grids_ = n_grids; 
        }
        inline void setRateConstants(std::vector<double>& rate_constants, size_t n_reactions){
            rate_constants_ = rate_constants.data(); 
            n_reactions_ = n_reactions; 
        }; 
        inline void setStateVariables(std::vector<double>& state_variables, size_t n_species){
            state_variables_ = state_variables.data(); 
            n_species_ = n_species; 
        }; 

        inline void setForcing(std::vector<double>& forcing, size_t n_species){
            forcing_ = forcing.data(); 
            n_species_ = n_species; 
        }; 
        
        inline void setJacobian(std::vector<double>& jacobian, size_t jacobian_size){
            jacobian_ = jacobian.data(); 
            jacobian_size_ = jacobian_size; 
        }; 

    }; //end class
}//end micm
