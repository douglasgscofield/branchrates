#ifndef TRAITMATRIX_H
#define TRAITMATRIX_H

// #include "branchrates.h"
#include "TaxonMatrix.h"
#include <vector>
#include <iostream>
#include <cassert>

class TraitMatrix {
    public:
        static const int    max_num_states;
        static const int    notset = -1;

        typedef int         trait_type;
        typedef std::vector<trait_type>
                            trait_vector_type;
        typedef std::vector<trait_vector_type>
                            trait_matrix_type;
        
        typedef int         freq_type;

    private:
        TaxonMatrix::taxon_id_vector_type
                            taxon_ids;  // across matrix.vals[]
        static const trait_type
                            valid_states[];
        trait_vector_type   trait_valid_states;
        TaxonMatrix*        taxa;
        
        void                realign_trait_matrix(
                                const TaxonMatrix::taxon_id_vector_type& newtid);
        void                realign_trait_matrix_to_leaf_ids();
        void                transpose_trait_matrix(trait_matrix_type& traitmat);


    public:
        class TraitMatrixEntry {
            private:
                int         id;
                freq_type   freq;

            public:
                trait_vector_type
                            vals;

                int         get_id() const;
                void        set_id(int i);
                const trait_vector_type&
                            get_vals() const;
                void        set_vals(const trait_vector_type& v);
                void        init_vals(const trait_vector_type& v);
                trait_type  get_trait(const TraitMatrix& tm,
                                      TaxonMatrix::taxon_id_type tid);
                int         get_freq() const;
                void        set_freq(int f);
                void        increment_freq();
                void        print(std::ostream& os = std::cout) const;
        };

        std::vector<TraitMatrixEntry>
                            matrix;
        const trait_matrix_type&
                            get_trait_matrix() const;
        void                set_trait_matrix(const trait_matrix_type& traitmat);

        TraitMatrix();
        ~TraitMatrix();
        
        void              read_file(const std::string& filename);
        void              write_file(const std::string& filename) const;

        // include ability to initialize from file, as below
        // include ability to simulate data, based on a phylogenetic tree??
        // TraitMatrix(string filename, bool file_is_summary = true);
        int               map_taxon_id_to_col(TaxonMatrix::taxon_id_type tid) const;
        void              set_taxon_ids(const TaxonMatrix::taxon_id_vector_type& tidm);
        TaxonMatrix*      get_taxa() const;
        void              set_taxa(TaxonMatrix& tm);
        const TaxonMatrix::taxon_id_vector_type&
                          get_taxon_ids() const;

        void     fill(const trait_matrix_type& traitmat,
                      const std::vector<freq_type>& freqs);

        void     fill(const trait_matrix_type& traitmat,
                      const std::vector<freq_type>& freqs,
                      const std::vector<std::string>& names);

        const trait_vector_type& get_trait_valid_states() const;
        
        void print(std::ostream& os = std::cout) const;
        friend std::ostream& operator<<(std::ostream& os,
                                        const TraitMatrix& tm);

};

//  TraitMatrix traits;
//  // ...code that fills traits...
//  for (int i = 0; i < traits.matrix.size(); ++i) {
//      TraitMatrix::trait_vector_type tv = traits.matrix.get_vals();
//      for (int tt = 0; tt < tv.size(); ++tt) {
//          // ... do something with tv[tt] ...
//      }
//      // or, do get the trait value for a particular taxa:
//      trait_type trait = traits.matrix[i].get_trait(traits, taxon_id);
//  }

///////////////////////////////////
/////////////////////////////////// TraitMatrix methods
///////////////////////////////////

inline TraitMatrix::TraitMatrix() : taxa(NULL) {
    /* empty */
}

inline TraitMatrix::~TraitMatrix() {
    /* empty */
}

inline const TaxonMatrix::taxon_id_vector_type& TraitMatrix::get_taxon_ids() const {
    return(taxon_ids);
}

inline void TraitMatrix::set_taxon_ids(
        const TaxonMatrix::taxon_id_vector_type& tidm) {
    taxon_ids = tidm;
}

inline TaxonMatrix* TraitMatrix::get_taxa() const {
    return(taxa);
}

inline const TraitMatrix::trait_vector_type& TraitMatrix::get_trait_valid_states() 
                                                                    const {
    return(trait_valid_states);
}

///////////////////////////////////
/////////////////////////////////// TraitMatrix::TraitMatrixEntry methods
///////////////////////////////////

inline TraitMatrix::trait_type TraitMatrix::TraitMatrixEntry::get_trait(
        const TraitMatrix& tm, TaxonMatrix::taxon_id_type tid) {
    int col = tm.map_taxon_id_to_col(tid);
    return(vals[col]);
}

inline int TraitMatrix::TraitMatrixEntry::get_id() const {
    return(id);
}

inline void TraitMatrix::TraitMatrixEntry::set_id(int i) {
    id = i;
}

inline const TraitMatrix::trait_vector_type&
             TraitMatrix::TraitMatrixEntry::get_vals() const {
    return(vals);
}

inline void TraitMatrix::TraitMatrixEntry::set_vals(
            const TraitMatrix::trait_vector_type& v) {
    assert(v.size() == vals.size());
    init_vals(v);
}

inline void TraitMatrix::TraitMatrixEntry::init_vals(
       const TraitMatrix::trait_vector_type& v) {
    vals = v;
}

inline TraitMatrix::freq_type TraitMatrix::TraitMatrixEntry::get_freq() const {
    return(freq);
}

inline void TraitMatrix::TraitMatrixEntry::increment_freq() {
    ++freq;
}

inline void TraitMatrix::TraitMatrixEntry::set_freq(TraitMatrix::freq_type i) {
    freq = i;
}

#endif // TRAITMATRIX_H

