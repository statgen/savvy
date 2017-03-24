#ifndef LIBVC_CMF_INDEX_READER_HPP
#define LIBVC_CMF_INDEX_READER_HPP

#include "portable_endian.hpp"

#include <fstream>
#include <cstring>
#include <stack>
#include <vector>
#include <algorithm>
#include <iostream>
#include <errno.h>

namespace vc
{
  namespace s1r
  {
    enum class block_size : std::uint16_t
    {
      bs_64 = 64,
      bs_512 = 512,
      bs_4096 = 4096,
      bs_32768 = 32768
    };

    class index_base
    {
    public:
      class internal_entry
      {
      public:
        std::uint64_t region_start() const { return be64toh(region_start_); }
        std::uint64_t region_length() const { return be64toh(region_length_); }
        std::uint64_t region_end() const { return this->region_start() + this->region_length(); }
        internal_entry() :
          region_start_(0),
          region_length_(0)
        { }
        internal_entry(std::uint64_t beg, std::uint64_t end) :
          region_start_(htobe64(beg)),
          region_length_(end > beg ? htobe64(end - beg) : 0)
        { }
      private:
        std::uint64_t region_start_;
        std::uint64_t region_length_;
      };

      class entry : public internal_entry
      {
      public:
        std::pair<std::uint64_t, std::uint64_t> value() const { return std::make_pair(be64toh(value1_), be64toh(value2_)); }

        entry() :
          value1_(0),
          value2_(0)
        { }

        entry(std::uint64_t beg, std::uint64_t end, std::uint64_t val1, std::uint64_t val2) :
          internal_entry(beg, end),
          value1_(htobe64(val1)),
          value2_(htobe64(val2))
        { }
      private:
        std::uint64_t value1_;
        std::uint64_t value2_;
      };

      struct node_position
      {
        node_position(std::uint64_t lev, std::uint64_t node_off) : level(lev), node_offset(node_off) {}
        std::uint64_t level;
        std::uint64_t node_offset;
      };

      struct tree_position : public node_position
      {
        tree_position(std::uint64_t lev, std::uint64_t node_off, std::uint64_t entry_off) : node_position(lev, node_off), entry_offset(entry_off) {}
        tree_position(node_position node_pos, std::uint64_t entry_off) : node_position(node_pos), entry_offset(entry_off) {}
        bool operator==(const tree_position& other) const { return level == other.level && node_offset == other.node_offset && entry_offset == other.entry_offset; }
        bool operator!=(const tree_position& other) const { return !(*this == other); }
        std::uint64_t entry_offset;
      };

      std::uint64_t header_block_count() const { return block_count_needed_for_header_; }
      std::uint64_t entry_count() const { return entry_count_; }
      std::uint64_t bucket_size() const { return block_size_; }
      std::uint64_t entry_count_at_level(std::size_t level) const { return entry_counts_per_level_[level]; }
      std::uint16_t entries_per_leaf_node() const { return block_size_ / std::uint16_t(32); }
      std::uint16_t entries_per_internal_node() const { return block_size_ / std::uint16_t(16); }

      std::uint64_t calculate_node_size(node_position input) const
      {
        std::uint64_t ret = 0;

        if (input.level < entry_counts_per_level_.size())
        {
          std::uint64_t entries_per_node = (input.level + 1 ==  entry_counts_per_level_.size() ? this->entries_per_leaf_node() : this->entries_per_internal_node());

          if (input.node_offset < entry_counts_per_level_[input.level] / entries_per_node)
            ret = entries_per_node;
          else
            ret = entry_counts_per_level_[input.level] % entries_per_node;
        }

        return ret;
      }

      tree_position calculate_parent_position(node_position input) const
      {
        return tree_position(input.level - 1, input.node_offset / this->entries_per_internal_node(), input.node_offset % this->entries_per_internal_node());
      }

      node_position calculate_child_position(tree_position input) const
      {
        return node_position(input.level + 1, this->entries_per_internal_node() * input.node_offset + input.entry_offset);
      }

      std::streampos calculate_file_position(node_position input) const
      {
        std::uint64_t ret = block_count_needed_for_header_ * block_size_;

        if (input.level > 0)
        {
          ret += block_size_; // for root node

          for (auto it = entry_counts_per_level_.begin(); it != entry_counts_per_level_.end(); ++it)
          {
            if (1 + std::distance(entry_counts_per_level_.begin(), it) < input.level)
            {
              ret += *it * block_size_;
            }
            else
            {
              ret += input.node_offset * block_size_;
              break;
            }
          }
        }

        return std::streampos(ret);
      }


      std::uint64_t tree_height() const
      {
        return entry_counts_per_level_.size();
      }
    protected:
      std::uint16_t block_count_needed_for_header_;
      std::uint16_t block_size_;
      std::vector<std::uint64_t> entry_counts_per_level_;
      std::uint64_t entry_count_;

      void init()
      {
        const std::size_t header_data_size = 7 + 2 + 8;
        block_count_needed_for_header_ = ceil_divide(header_data_size, block_size_);

        this->entry_counts_per_level_.insert(this->entry_counts_per_level_.begin(), entry_count_);
        for (std::uint64_t nodes_at_current_level = ceil_divide(entry_count_, (std::uint64_t) this->entries_per_leaf_node());
             nodes_at_current_level > 1;
             nodes_at_current_level = ceil_divide(nodes_at_current_level, (std::uint64_t) this->entries_per_internal_node()))
        {
          this->entry_counts_per_level_.insert(this->entry_counts_per_level_.begin(), nodes_at_current_level);
        }
      }

      static std::uint64_t ceil_divide(std::uint64_t x, std::uint64_t y)
      {
        return (x + y - 1) / y;
      }

      static std::uint16_t ceil_divide(std::uint16_t x, std::uint16_t y)
      {
        return (x + y - std::uint16_t(1)) / y;
      }
    };

    class reader : public index_base
    {
    public:
      class query
      {
      public:
        class iterator
        {
        public:
          typedef iterator self_type;
          typedef std::ptrdiff_t difference_type;
          typedef entry value_type;
          typedef const value_type& reference;
          typedef const value_type* pointer;
          typedef std::bidirectional_iterator_tag iterator_category;

          iterator(reader& rdr, std::istream& ifs, std::uint64_t beg, std::uint64_t end, tree_position pos = {0, 0, 0}) :
            reader_(&rdr),
            ifs_(&ifs),
            beg_(beg),
            end_(end),
            leaf_node_(reader_->entries_per_leaf_node()),
            position_(pos)
          {
            if (position_ != end_tree_position())
            {
              const std::uint64_t leaf_level = reader_->tree_height() - 1;
              ifs_->seekg(reader_->calculate_file_position(position_));
              if (ifs_->eof())
              {
                std::cout << ifs_->good() << std::endl;
              }
              if (position_.level == leaf_level)
                ifs_->read((char*) leaf_node_.data(), reader_->block_size_);
              else
              {
                traversal_chain_.emplace(reader_->entries_per_internal_node());
                ifs_->read((char*) (traversal_chain_.top().data()), reader_->block_size_);
              }

              traverse_right();
            }
          }

          void traverse_right()
          {
            const std::uint64_t leaf_level = reader_->tree_height() - 1;
            bool leaf_entry_found = false;
            while (!leaf_entry_found && position_ != end_tree_position())
            {
              if (position_.level == leaf_level)
              {
                auto entry_it = leaf_node_.begin() + position_.entry_offset;
                const auto entry_end_it = leaf_node_.begin() + reader_->calculate_node_size(position_);
                for (; entry_it != entry_end_it; ++entry_it)
                {
                  if (entry_it->region_start() <= end_ && entry_it->region_end() >= beg_)
                    break;
                }

                if (entry_it == entry_end_it)
                {
                  position_ = reader_->calculate_parent_position(position_);
                  position_ = tree_position(position_, position_.entry_offset + 1);
                }
                else
                {
                  position_ = tree_position(position_.level, position_.node_offset, (std::uint64_t) std::distance(leaf_node_.begin(), entry_it));
                  leaf_entry_found = true;
                }
              }
              else
              {
                auto entry_it = traversal_chain_.top().begin() + position_.entry_offset;
                const auto entry_end_it = traversal_chain_.top().begin() + reader_->calculate_node_size(position_);
                for (; entry_it != entry_end_it; ++entry_it)
                {
                  if (entry_it->region_start() <= end_ && entry_it->region_end() >= beg_)
                    break;
                }

                if (entry_it == entry_end_it)
                {
                  position_ = reader_->calculate_parent_position(position_);
                  position_ = tree_position(position_, position_.entry_offset + 1);
                  traversal_chain_.pop();
                }
                else
                {
                  position_ = tree_position(reader_->calculate_child_position(tree_position(position_, (std::uint64_t) std::distance(traversal_chain_.top().begin(), entry_it))), 0);

                  ifs_->seekg(reader_->calculate_file_position(position_));
                  if (position_.level == leaf_level)
                    ifs_->read((char*) leaf_node_.data(), reader_->block_size_);
                  else
                  {
                    traversal_chain_.emplace(reader_->entries_per_internal_node());
                    ifs_->read((char*) (traversal_chain_.top().data()), reader_->block_size_);
                  }
                }
              }
            }
          }

          self_type& operator++()
          {
            position_ = tree_position(position_, position_.entry_offset + 1);
            traverse_right();
            return *this;
          }

          self_type operator++(int)
          {
            self_type r = *this;
            ++(*this);
            return r;
          }

          reference operator*() { return leaf_node_[position_.entry_offset]; }
          pointer operator->() { return &(leaf_node_[position_.entry_offset]); }
          bool operator==(const self_type& other) { return position_ == other.position_; }
          bool operator!=(const self_type& other) { return position_ != other.position_; }

        private:
          reader* reader_;
          std::istream* ifs_;
          std::uint64_t beg_;
          std::uint64_t end_;
          std::stack<std::vector<internal_entry>> traversal_chain_;
          std::vector<entry> leaf_node_;
          tree_position position_;

          static const tree_position end_tree_position() { return tree_position((std::uint64_t)-1, 0, 1); }
        };

        query(reader& rdr, std::istream& ifs, std::uint64_t beg, std::uint64_t end) :
          reader_(&rdr),
          ifs_(&ifs),
          beg_(beg),
          end_(end)
        {
        }

        iterator begin() { return iterator(*reader_, *ifs_, beg_, end_); }
        iterator end() { return iterator(*reader_, *ifs_, beg_, end_, tree_position((std::uint64_t)(-1), 0, 1)); }
      private:
        reader* reader_;
        std::istream* ifs_;
        std::uint64_t beg_;
        std::uint64_t end_;
      };

      reader(const std::string& file_path) :
        ifs_(file_path, std::ios::binary)
      {
        std::string version(7, '\0');
        ifs_.read(&version[0], version.size());
        std::uint8_t block_size_exponent;
        ifs_.read((char*)(&block_size_exponent), 1);
        block_size_ = static_cast<std::uint16_t>(std::pow(8.0, (0x03 & block_size_exponent) + 2));
        ifs_.read((char*)(&entry_count_), 8);
        entry_count_ = be64toh(entry_count_);

        index_base::init();
      }

      query create_query(std::uint64_t beg, std::uint64_t end)
      {
        query ret(*this, ifs_, beg, end);
        return ret;
      }
    private:
      std::ifstream ifs_;
    };

    template<typename EntryIter>
    class writer : public index_base
    {
    public:
      writer(const std::string& file_path, EntryIter beg, EntryIter end, block_size bs = block_size::bs_4096) :
        ofs_(file_path, std::ios::binary)
      {
        std::sort(beg, end, [](const reader::entry& a, const reader::entry& b)
        {
          double mid_a = static_cast<double>(a.region_start()) + (static_cast<double>(a.region_length()) / 2.0);
          double mid_b = static_cast<double>(b.region_start()) + (static_cast<double>(b.region_length()) / 2.0);

          return mid_a < mid_b;
        });

        std::uint8_t block_size_exponent;
        switch (bs)
        {
          case block_size::bs_64:
            block_size_ = 64;
            block_size_exponent = 2;
            break;
          case block_size::bs_512:
            block_size_ = 512;
            block_size_exponent = 3;
            break;
          case block_size::bs_4096:
            block_size_ = 4096;
            block_size_exponent = 4;
            break;
          case block_size::bs_32768:
            block_size_ = 32768;
            block_size_exponent = 5;
            break;
        }
        block_size_exponent -= 2;

        entry_count_ = (std::uint64_t) std::distance(beg, end);
        std::uint64_t entry_count_be = htobe64(entry_count_);

        index_base::init();

        std::string header_block(16, '\0');
        std::memcpy(&header_block[0], "s1r\0x00\0x01\0x00\0x00", 7);
        std::memcpy(&header_block[7], &block_size_exponent, 1);
        std::memcpy(&header_block[8], &entry_count_be, 8);

        header_block.resize(static_cast<std::uint16_t>(bs), '\0');
        ofs_.write(header_block.data(), header_block.size());

        std::vector<std::pair<std::vector<internal_entry>, tree_position>> current_nodes_at_each_internal_level;
        current_nodes_at_each_internal_level.reserve(this->tree_height() - 1);
        for (std::size_t i = 0; i < this->tree_height() - 1; ++i)
          current_nodes_at_each_internal_level.emplace_back(std::make_pair(std::vector<internal_entry>(this->entries_per_internal_node()), tree_position(i, 0, 0)));

        std::vector<entry> current_leaf_node(this->entries_per_leaf_node());
        tree_position current_leaf_position(this->tree_height() - 1, 0, 0);

        std::size_t node_counter = 0;
        for (auto entry_it = beg; entry_it != end; ++entry_it)
        {
          const bool last_entry = entry_it + 1 == end;
          reader::entry& e = *entry_it;

          current_leaf_node[current_leaf_position.entry_offset] = e;

          if (current_leaf_position.entry_offset + 1 == current_leaf_node.size() || last_entry)
          {
            ofs_.seekp(this->calculate_file_position(current_leaf_position));
            ofs_.write((char*)(current_leaf_node.data()), block_size_);
            ++node_counter;

            std::uint64_t node_range_min = (std::uint64_t)-1;
            std::uint64_t node_range_max = 0;
            for (std::size_t i = 0; i <= current_leaf_position.entry_offset; ++i)
            {
              node_range_min = std::min(node_range_min, current_leaf_node[i].region_start());
              node_range_max = std::max(node_range_max, current_leaf_node[i].region_end());
            }

            for (auto rit = current_nodes_at_each_internal_level.rbegin(); rit != current_nodes_at_each_internal_level.rend(); ++rit)
            {
              rit->first[rit->second.entry_offset] = internal_entry(node_range_min, node_range_max);



              if (rit->second.entry_offset + 1 == rit->first.size() || last_entry)
              {
                ofs_.seekp(this->calculate_file_position(rit->second));
                ofs_.write((char*)(rit->first.data()), block_size_);
                ++node_counter;

                node_range_min = (std::uint64_t)-1;
                node_range_max = 0;
                for (std::size_t i = 0; i <= rit->second.entry_offset; ++i)
                {
                  node_range_min = std::min(node_range_min, rit->first[i].region_start());
                  node_range_max = std::max(node_range_max, rit->first[i].region_end());
                }

                rit->first = std::vector<internal_entry>(this->entries_per_internal_node());
                ++(rit->second.node_offset);
                rit->second.entry_offset = 0;
              }
              else
              {
                ++(rit->second.entry_offset);
                break;
              }
            }

            current_leaf_node = std::vector<entry>(this->entries_per_leaf_node());
            ++(current_leaf_position.node_offset);
            current_leaf_position.entry_offset = 0;
          }
          else
          {
            ++(current_leaf_position.entry_offset);
          }
        }
      }

      bool flush()
      {
        return ofs_.flush().good();
      }
    private:
      std::ofstream ofs_;
    };

    template<typename EntryIter>
    bool create_file(const std::string& file_path, EntryIter beg, EntryIter end, block_size bs = block_size::bs_4096)
    {
      writer<EntryIter> w(file_path, beg, end, bs);
      return w.flush();
    }
  }
}

#endif //LIBVC_CMF_INDEX_READER_HPP