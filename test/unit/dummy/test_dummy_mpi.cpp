#include <micm/dummy.hpp>
#include <assert.h>
#include <mpi.h>
#include <vector>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  std::vector<bool> failed{false};

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  try{
    double d = dummy(1.1, my_rank);
    decl();
  }
  catch(...)
  {
    failed[my_rank] = true;
  }

  for(uint8_t i = 0; i < world_size; ++i)
  {
    assert(!failed[i]);
  }

  MPI_Finalize();
}