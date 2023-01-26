######################################################################
# MUSICA library

set(MUSICA_SRC component.F90
               component_set.F90
               config.F90
               domain.F90
               domain_factory.F90
               domain_iterator.F90
               domain_state.F90
               domain_state_accessor.F90
               domain_state_mutator.F90
               domain_target_cells.F90
               domain_target_columns.F90
               domain_target_model_top_cells.F90
               domain_target_surface_cells.F90
               domains/cell.F90
               evolving_conditions.F90
               grid.F90
               initial_conditions.F90
               input_output/file.F90
               input_output/file_dimension.F90
               input_output/file_dimension_range.F90
               input_output/file_dimension_factory.F90
               input_output/file_dimensions/netcdf.F90
               input_output/file_dimensions/text.F90
               input_output/file_factory.F90
               input_output/file_paired_variable.F90
               input_output/file_updater.F90
               input_output/file_variable.F90
               input_output/file_variable_factory.F90
               input_output/file_variables/netcdf.F90
               input_output/file_variables/text.F90
               input_output/files/netcdf.F90
               input_output/files/text.F90
               input_output/input_output_processor.F90
               interpolator.F90
               interpolator_strategies/linear_1D.F90
               io.F90
               io/netcdf.F90
               lookup/axis.F90
               lookup/2D_axis.F90
               util/array.F90
               util/assert.F90
               util/constants.F90
               util/convert.F90
               util/data_type.F90
               util/datetime.F90
               util/iterator.F90
               util/logger.F90
               util/math.F90
               util/property.F90
               util/property_set.F90
               util/string.F90
               util/target.F90
               )

add_library(musicalib ${MUSICA_SRC})

target_link_libraries(musicalib ${JSON_LIB} ${NETCDF_LIBS})

set_target_properties(musicalib PROPERTIES OUTPUT_NAME musica)

######################################################################
