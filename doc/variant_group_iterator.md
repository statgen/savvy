# Variant Group Iterator

## Collapsed Dosage
```c++
#include <savvy/reader.hpp>
#include <savvy/variant_group_iterator.hpp>

std::ifstream marker_group_file("file_groups.txt", std::ios::binary);
savvy::indexed_reader marker_file("marker_file.sav", {""}, savvy::fmt::ds);
std::string marker_group_line;
std::vector<float> collapsed_dose(marker_file.samples().size());

while (std::getline(marker_group_file, marker_group_line))
{
  savvy::variant_group_iterator<savvy::compressed_vector<float>> it(marker_file, marker_group_line);
  savvy::variant_group_iterator<savvy::compressed_vector<float>> end{};

  std::fill(collapsed_dose.begin(), collapsed_dose.end(), 0.f);

  for ( ; it != end; ++it)
  {
    it.group_id();
    it.sites();

    std::string marker_id = it->chromosome() + ":" + std::to_string(it->position()) + "_" + it->ref() + "/" + it->alt();
    
    for (auto dose_it = it->data().begin(); dose_it != it->data().end(); ++dose_it)
      collapsed_dose[dose_it.offset()] += *dose_it;
  }
}
```

## Variant-Sample Matrix
```c++
#include <savvy/reader.hpp>
#include <savvy/variant_group_iterator.hpp>

std::ifstream marker_group_file("file_groups.txt", std::ios::binary);
savvy::indexed_reader marker_file("marker_file.sav", {""}, savvy::fmt::ds);
std::string marker_group_line;
std::size_t sample_size = marker_file.samples().size();
std::vector<float> group_matrix;

while (std::getline(marker_group_file, marker_group_line))
{
  savvy::variant_group_iterator<savvy::compressed_vector<float>> it(marker_file, marker_group_line);
  savvy::variant_group_iterator<savvy::compressed_vector<float>> end{};

  group_matrix.resize(0);
  if (it != end)
    group_matrix.resize(it.sites().size() * sample_size);

  std::size_t cnt = 0;
  for ( ; it != end; ++it,++cnt)
  {
    it.group_id();
    it.sites();

    std::string marker_id = it->chromosome() + ":" + std::to_string(it->position()) + "_" + it->ref() + "/" + it->alt();
    
    for (auto dose_it = it->data().begin(); dose_it != it->data().end(); ++dose_it)
      group_matrix[cnt * sample_size + dose_it.offset()] = *dose_it;
  }
  
  group_matrix.resize(sample_size * cnt);
}
```
