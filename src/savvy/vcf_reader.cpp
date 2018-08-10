/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "savvy/vcf_reader.hpp"
extern "C" {
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>
}
#include <shrinkwrap/gz.hpp>

#include <sstream>
#include <limits>

namespace savvy
{
  namespace vcf
  {
    // https://github.com/samtools/htslib/blob/master/tabix.c#L229-L262
    std::vector<std::string> query_chromosomes(const std::string& file_path)
    {
      std::vector<std::string> ret;

      if (::savvy::detail::has_extension(file_path, ".vcf") || ::savvy::detail::has_extension(file_path, ".vcf.gz"))
      {
        tbx_t *tbx = tbx_index_load(file_path.c_str());
        if (tbx)
        {
          const char **seq;
          int nseq{};
          seq = tbx_seqnames(tbx, &nseq);
          ret.resize(nseq);
          for (int i = 0; i < nseq; ++i)
            ret[i] = seq[i];
          tbx_destroy(tbx);
          free(seq);
        }
      }
      else if (::savvy::detail::has_extension(file_path, ".bcf"))
      {
        htsFile *fp = hts_open(file_path.c_str(),"r");
        if (fp)
        {
          bcf_hdr_t *hdr = bcf_hdr_read(fp);
          if (hdr)
          {
            hts_idx_t *idx = bcf_index_load(file_path.c_str());
            if (idx)
            {
              const char **seq;
              int nseq{};
              seq = bcf_index_seqnames(idx, hdr, &nseq);
              ret.resize(nseq);
              for (int i = 0; i < nseq; ++i)
                ret[i] = seq[i];
              hts_idx_destroy(idx);
              free(seq);
            }
            bcf_hdr_destroy(hdr);
          }
          hts_close(fp);
        }
      }

      return ret;
    }

    class hts_file : public detail::hts_file_base
    {
    public:
      hts_file(htsFile* fp, bcf_hdr_t* hdr, bcf1_t* rec, bool destroy_hdr = true) :
        file_(fp),
        hdr_(hdr),
        rec_(rec),
        destroy_hdr_(destroy_hdr)
      {
      }

      virtual ~hts_file()
      {
        try
        {
          if (hdr_ && destroy_hdr_)
            bcf_hdr_destroy(hdr_);
          if (file_)
            bcf_close(file_);
          if (rec_ && destroy_hdr_)
            bcf_destroy1(rec_);
        }
        catch (...)
        {
          assert(!"bcf_hdr_destroy or bcf_close threw exception!");
        }
      }

      void init_headers(std::vector<std::pair<std::string,std::string>>& destination)
      {
        if (hdr_)
        {
          destination.reserve(std::size_t(hdr_->nhrec - 1));
          for (int i = 0; i < hdr_->nhrec; ++i)
          {
            std::string key, val;
            if (hdr_->hrec[i]->key && hdr_->hrec[i]->value)
            {
              key = hdr_->hrec[i]->key;
              val = hdr_->hrec[i]->value;
            }
            else if (hdr_->hrec[i]->key && hdr_->hrec[i]->nkeys) // (hdr->hrec[i]->type == BCF_HL_INFO || hdr->hrec[i]->type == BCF_HL_FLT || hdr->hrec[i]->type == BCF_HL_STR))
            {
              bcf_hrec_t* r = hdr_->hrec[i];
              key = r->key;
              std::stringstream ss_val;

              ss_val << "<";
              for (int j = 0; j < r->nkeys - 1; ++j) // minus 1 to remove IDX;
              {
                if (j > 0)
                  ss_val << ",";
                if (r->keys[j])
                  ss_val << r->keys[j];
                ss_val << "=";
                if (r->vals[j])
                  ss_val << r->vals[j];
              }
              ss_val << ">";
              val = ss_val.str();
            }

            if (key.size())
              destination.emplace_back(std::move(key), std::move(val));
            //ret.insert(std::upper_bound(ret.begin(), ret.end(), std::make_pair(key, std::string()), [](const auto& a, const auto& b) { return a.first < b.first; }), {std::move(key), std::move(val)});
          }
        }
      }

      void init_sample_ids(std::vector<std::string>& destination)
      {
        if (hdr_)
        {
          const char **beg = (const char **) (hdr_->samples);
          const char **end = (const char **) (hdr_->samples + bcf_hdr_nsamples(hdr_));
          destination.reserve(end - beg);

          for (; beg != end; ++beg)
          {
            destination.emplace_back(*beg);
          }
        }
      }

      void init_info_fields(std::vector<std::string>& destination)
      {
        if (hdr_)
        {
          destination = {"ID", "QUAL", "FILTER"};
          std::unordered_set<std::string> unique_info_fields{destination.begin(), destination.end()};
          for (int i = 0; i < hdr_->nhrec; ++i)
          {
            if (hdr_->hrec[i]->type == BCF_HL_INFO)
            {
              bcf_hrec_t* r = hdr_->hrec[i];
              for (int j = 0; j < r->nkeys; ++j)
              {
                if (strcmp(r->keys[j], "ID") == 0)
                {
                  const char* inf = r->vals[j];
                  if (inf && unique_info_fields.emplace(inf).second)
                    destination.emplace_back(inf);
                }
              }
            }
          }
        }
      }

      std::size_t cur_fmt_field_size() const
      {
        if (rec_)
          return rec_->n_fmt;
        return 0;
      }

      const char*const cur_fmt_field(std::size_t idx) const
      {
        int fmt_id = rec_->d.fmt[idx].id;
        return hdr_->id[BCF_DT_ID][fmt_id].key;
      }

      std::size_t cur_num_alleles() const
      {
        if (rec_)
          return rec_->n_allele;
        return 0;
      }

      virtual bool read_next_record()
      {
        if (bcf_read(file_, hdr_, rec_) >= 0)
        {
          bcf_unpack(rec_, BCF_UN_ALL);
          return true;
        }
        return false;
      }

      bool get_cur_format_values_int32(const char* tag, int**buf, int*sz) const
      {
        return bcf_get_format_int32(hdr_, rec_, tag, buf, sz) >= 0;
      }

      bool get_cur_format_values_float(const char* tag, int**buf, int*sz) const
      {
        return bcf_get_format_float(hdr_, rec_, tag, buf, sz) >= 0;
      }

      site_info cur_site_info(std::size_t allele_index) const
      {
        std::size_t n_info = rec_->n_info;
        std::size_t n_flt = rec_->d.n_flt;
        bcf_info_t* info = rec_->d.info;
        std::unordered_map<std::string, std::string> props;
        props.reserve(n_info + 2);

        if (std::isnan(rec_->qual))
        {
          props["QUAL"] = ".";
        }
        else
        {
          std::string qual(std::to_string(rec_->qual));
          qual.erase(qual.find_last_not_of('0') + 1); // rtrim zeros.
          qual.erase(qual.find_last_not_of('.') + 1);
          props["QUAL"] = std::move(qual);
        }

        std::stringstream ss;
        for (std::size_t i = 0; i < n_flt; ++i)
        {
          if (i > 0)
            ss << ";";
          ss << bcf_hdr_int2id(hdr_, BCF_DT_ID, rec_->d.flt[i]);
        }
        std::string fltr(ss.str());
        if (fltr == ".")
          fltr.clear();
        props["FILTER"] = std::move(fltr);
        props["ID"] = rec_->d.id;


        for (std::size_t i = 0; i < n_info; ++i)
        {
          // bcf_hdr_t::id[BCF_DT_ID][$key].key
          const char* key = hdr_->id[BCF_DT_ID][info[i].key].key;
          if (key)
          {
            switch (info[i].type)
            {
            case BCF_BT_NULL:
              props[key] = "1"; // Flag present so should be true.
              break;
            case BCF_BT_INT8:
            case BCF_BT_INT16:
            case BCF_BT_INT32:
              props[key] = std::to_string(info[i].v1.i);
              break;
            case BCF_BT_FLOAT:
              props[key] = std::to_string(info[i].v1.f);
              props[key].erase(props[key].find_last_not_of('0') + 1); // rtrim zeros.
              props[key].erase(props[key].find_last_not_of('.') + 1);
              break;
            case BCF_BT_CHAR:
              props[key] = std::string((char*)info[i].vptr, info[i].vptr_len);
              break;
            }
          }
        }

        return site_info(
          std::string(bcf_hdr_id2name(hdr_, rec_->rid)),
          static_cast<std::uint64_t>(rec_->pos + 1),
          std::string(rec_->d.allele[0]),
          std::string(rec_->n_allele > 1 ? rec_->d.allele[allele_index] : ""),
          std::move(props));
      }
    protected:
      bcf_hdr_t* hdr_ = nullptr;
      bcf1_t* rec_ = nullptr;
    private:
      htsFile* file_ = nullptr;
      bool destroy_hdr_ = false;
    };

    class hts_indexed_file : public hts_file
    {
    public:
      hts_indexed_file(bcf_srs_t* sr, bcf_hdr_t* hdr) :
        hts_file(nullptr, hdr, nullptr, false),
        synced_readers_(sr)
      {
      }

      ~hts_indexed_file()
      {
        if (synced_readers_)
          bcf_sr_destroy(synced_readers_);
      }

      bool read_next_record()
      {
        if (bcf_sr_next_line(synced_readers_) && (rec_ = bcf_sr_get_line(synced_readers_, 0)))
        {
          bcf_unpack(rec_, BCF_UN_ALL);
          return true;
        }
        return false;
      }
    private:
      bcf_srs_t* synced_readers_ = nullptr;
    };

    std::unique_ptr<detail::hts_file_base> detail::hts_file_base::create_file(const std::string& file_path)
    {
      htsFile* fp = bcf_open(file_path.c_str(), "r");
      bcf_hdr_t* hdr = nullptr;
      bcf1_t* rec = bcf_init1();

      if (fp && rec)
      {
        hdr = bcf_hdr_read(fp);
        if (hdr)
        {
          return std::unique_ptr<detail::hts_file_base>(new hts_file(fp, hdr, rec));
        }
      }


      if (fp)
        bcf_close(fp);
      if (rec)
        bcf_destroy1(rec);

      return nullptr;
    }

    std::unique_ptr<detail::hts_file_base> detail::hts_file_base::create_indexed_file(const std::string& file_path, const region& reg)
    {
      bcf_srs_t* sr = bcf_sr_init();
      //bcf1_t* rec = bcf_init1();

      if (sr)// && rec)
      {
        std::stringstream contigs;
        contigs << reg.chromosome();
        if (reg.from() > 1 || reg.to() != std::numeric_limits<std::uint64_t>::max())
          contigs << ":" << reg.from() << "-" << reg.to();

        if (bcf_sr_set_regions(sr, contigs.str().c_str(), 0) == 0 && bcf_sr_add_reader(sr, file_path.c_str()) == 1)
        {
          bcf_hdr_t* hdr = bcf_sr_get_header(sr, 0);
          if (hdr)
          {
            return std::unique_ptr<detail::hts_file_base>(new hts_indexed_file(sr, hdr));
          }
        }
      }

      if (sr)
        bcf_sr_destroy(sr);

//      if (rec)
//        bcf_destroy1(rec);

      return nullptr;
    }

    std::unique_ptr<std::ostream> detail::create_out_stream(const std::string& file_path, compression_type type)
    {
      if (type == compression_type::none)
        return std::unique_ptr<std::ostream>(new std::ofstream(file_path));
      else
        return std::unique_ptr<std::ostream>(new shrinkwrap::bgzf::ostream(file_path));
    }
  }
}

//savvy::vcf::marker::const_iterator operator+(savvy::vcf::marker::const_iterator::difference_type n, const savvy::vcf::marker::const_iterator& a)
//{
//  savvy::vcf::marker::const_iterator ret(a);
//  return (ret += n);
//}
