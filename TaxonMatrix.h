#ifndef TAXON_H
#define TAXON_H

#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>

//////////////////////
// Below these comments are data structures which define internal taxa
// that are used to establish a taxon matrix for testing purposes
//////////////////////
struct s_init_matrix {
    int           id;
    const char*   name;
    const char*   description;
};

static s_init_matrix init_matrix[] = {
    {0,  "a", ""},
    {0,  "b", ""},
    {0,  "c", ""},
    {-1, "",  ""},
};
//////////////////////
// Above these comments are data structures which define internal taxa
// that are used to establish a taxon matrix for testing purposes
//////////////////////

class TaxonMatrix {
    public:
        typedef int                          taxon_id_type;
        typedef std::vector<taxon_id_type>   taxon_id_vector_type;
        static const int                     notset;
        
        class Taxon {
            private:
                taxon_id_type      id;
                std::string        name;
                std::string        description;
                int                leaf_id;
                
            public:
                const std::string& get_name() const;
                void               set_name(const std::string& n);
                const std::string& get_description() const;
                void               set_description(const std::string& d);
                taxon_id_type      get_id() const;
                void               set_id(taxon_id_type i);
                int                get_leaf_id() const;
                void               set_leaf_id(int lid);
                void               print(std::ostream& os = std::cout) const;

        };
        typedef std::vector<Taxon> taxon_vector_type;
        
        taxon_vector_type          taxa;

        TaxonMatrix();
        TaxonMatrix(const std::string& filename);
        ~TaxonMatrix();

        void                       read_file(const std::string& filename);
        void                       read_internal_data();
        taxon_id_type              taxon_id_for_leaf_id(int lid) const;
        taxon_id_type              add_taxon(const std::string& n,
                                             const std::string& d = "");
        taxon_id_type              get_taxon_id(const std::string& name) const;

        void                       print(std::ostream& os = std::cout) const;

        friend std::ostream&       operator<<(std::ostream& os,
                                              const TaxonMatrix& tm);
};

/////////////////////////////////////////
///////////////////////////////////////// TaxonMatrix::Taxon methods
/////////////////////////////////////////

inline int TaxonMatrix::Taxon::get_leaf_id() const {
    return(leaf_id);
}

inline void TaxonMatrix::Taxon::set_leaf_id(int lid) {
    leaf_id = lid;
}

inline const std::string& TaxonMatrix::Taxon::get_name() const {
    return(name);
}

inline void TaxonMatrix::Taxon::set_name(const std::string& n) {
    name = n;
}

inline const std::string& TaxonMatrix::Taxon::get_description() const {
    return(description);
}

inline void TaxonMatrix::Taxon::set_description(const std::string& d) {
    description = d;
}

inline TaxonMatrix::taxon_id_type TaxonMatrix::Taxon::get_id() const {
    return(id);
}

inline void TaxonMatrix::Taxon::set_id(TaxonMatrix::taxon_id_type i) {
    id = i;
}

/////////////////////////////////////////
///////////////////////////////////////// TaxonMatrix methods
/////////////////////////////////////////

inline TaxonMatrix::TaxonMatrix() {
    /* empty */
}

inline TaxonMatrix::TaxonMatrix(const std::string& filename) {
    read_file(filename);
}

inline TaxonMatrix::~TaxonMatrix() {
    /* empty */
}

inline TaxonMatrix::taxon_id_type TaxonMatrix::taxon_id_for_leaf_id(int lid) const
{
    assert(lid >= 0);
    int i;
    for (i = 0; i < taxa.size(); ++i)
        if (taxa[i].get_leaf_id() == lid)
            break;
    if (i >= taxa.size())
        return(notset);
    else
        return(i);
}


inline void TaxonMatrix::read_internal_data()
{
    for (int i = 0; init_matrix[i].id >= 0; ++i) {
        // we don't use the id field, and probably never will
        // we also don't need to know the id that resulted from the add_taxon()
        TaxonMatrix::taxon_id_type id = add_taxon(init_matrix[i].name,
                                                  init_matrix[i].description);
    }
}

inline TaxonMatrix::taxon_id_type TaxonMatrix::add_taxon(const std::string& n,
                                                         const std::string& d) {
    for (int i = 0; i < taxa.size(); ++i) {
        assert(taxa[i].get_name() != n);
    }
    int end = taxa.size();
    taxa.resize(end + 1);
    taxa[end].set_id(end);
    taxa[end].set_name(n);
    taxa[end].set_description(d);
    taxa[end].set_leaf_id(notset);
    return(end);
}

inline TaxonMatrix::taxon_id_type TaxonMatrix::get_taxon_id(
                                                   const std::string& n) const
{
    for (int i = 0; i < taxa.size(); ++i) {
        if(taxa[i].get_name() == n) {
            return(taxa[i].get_id());
        }
    }
    return(notset);
}


#endif // TAXON_H

