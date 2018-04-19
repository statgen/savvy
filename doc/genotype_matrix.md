# Targeting Matrices

BCF, VCF and SAV file formats all support having variants of different ploidy levels in the same file. This freedom makes it difficult to have built-in support for statically sized 3rd-party matrix classes in the Savvy library since some data formats (e.g., allele, haplotype_dosage, etc.) depend on ploidy level. Applications, on the other hand, can make consistent ploidy level a restriction, and matrix integration with Savvy can be done with a small wrapper class. 

The wrapper class needs to define a resize method and a subscript operator. Every call to read methods in Savvy will resize the the target vector to zero before resizing it to the correct size for next variant to be read. The following example shows how to do this with an Armadillo sparse matrix, but the same technique can be used for virtually any matrix class.

   

```c++
#include <armadillo>
#include <system_error>
#include <savvy/reader.hpp>

template <typename T>
class sparse_matrix_wrapper
{
public:
  typedef T value_type;

  sparse_matrix_wrapper(arma::SpMat<T>& target_matrix, std::size_t target_row, std::error_code& ec) :
    mat_(target_matrix),
    ec_(ec),
    row_index_(target_row)
  {}

  void resize(std::size_t sz)
  {
    if (sz == 0)
    {
      // Do nothing.
    }
    else if (sz == mat_.n_cols)
    {
      // Do nothing.
    }
    else
    {
      ec_ = std::make_error_code(std::errc::argument_out_of_domain);
    }
  }

  arma::MapMat_elem<T> operator[](std::size_t idx)
  {
    return mat_(row_index_, idx);
  }

  sparse_matrix_wrapper& ref() { return *this; }
private:
  arma::SpMat<T>& mat_;
  std::error_code& ec_;
  std::size_t row_index_;
};
```

This wrapper class can then be used as seen below.

```c++
int main()
{
  const std::size_t window_size = 100;
  savvy::reader input_file("file.sav", savvy::fmt::gt);
  std::vector<savvy::site_info> annotations(window_size);
  arma::SpMat<float> genotype_matrix(window_size, input_file.samples().size() * 2);
  
  std::error_code ec;
  std::size_t row = 0;
  while (input_file.read(annotations[row], sparse_matrix_wrapper<float>(genotype_matrix, row, ec).ref()))
  {
    if (ec)
    {
      // ploidy not 2 / handle error
    }
    else
    {
      ++row;
      if (row == genotype_matrix.n_rows)
      {
        // process full genotype matrix ...
  
        genotype_matrix.zeros(); // reset matrix to all zeros.
        row = 0; // reset row counter.
      }
    }
  }
  
  if (row)
  {
    // process partial genotype matrix ...
  }
  
  return 0;
}
```

If your matrix of choice is a vector of vectors, then this can be done without a wrapper class.

```c++
int main()
{
  savvy::reader input_file("/Users/lefaivej/Developer/projects/savvy/test_file.vcf", savvy::fmt::gt);
  const std::size_t window_size = 100;
  const std::size_t num_columns = input_file.samples().size() * 2;
  std::vector<savvy::site_info> annotations(window_size);
  std::vector<std::vector<float>> genotype_matrix(window_size);

  std::size_t row = 0;
  while (input_file.read(annotations[row], genotype_matrix[row]))
  {
    if (genotype_matrix[row].size() != num_columns)
    {
      // ploidy not 2 / handle error
    }
    else
    {
      ++row;
      if (row == genotype_matrix.size())
      {
        // process full genotype matrix ...

        row = 0; // reset row counter
      }
    }
  }

  if (row)
  {
    // process partial genotype matrix ...
  }
  
  return 0;
}
``` 