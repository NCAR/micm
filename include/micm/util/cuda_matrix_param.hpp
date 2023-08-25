#pragma once
#include <vector> 
namespace micm{
    class CUDAMatrixParam{
        public: 
        const double* rate_constants_; 
        const double* state_variables_; 
        double* forcing_; 
        double* jacobian_; 
        const size_t n_grids_; 
        const size_t n_reactions_; 
        const size_t n_species_; 
        size_t jacobian_size_; 
        bool set_n_species ; 
        
        CUDAMatrixParam(){
            set_n_species = false; 
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
            if (set_n_species == false){
            n_species_ = n_species; 
            set_n_species = true; 
            }
        }; 
        inline void setForcing(std::vector<double>& forcing, size_t n_species){
            forcing_ = forcing.data(); 
            if (set_n_species == false){
            n_species_ = n_species; 
            set_n_species = true; 
            }
        }; 
        inline void setJacobian(std::vector<double>& jacobian, size_t jacobian_size){
            jacobian_ = jacobian.data(); 
            jacobian_size_ = jacobian_size; 
        }

    }; //end class
}//end micm
