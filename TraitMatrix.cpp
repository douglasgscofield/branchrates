#include "TraitMatrix.h"

const int                       TraitMatrix::max_num_states = 2;
const TraitMatrix::trait_type   TraitMatrix::valid_states[] = {0, 1};


const TraitMatrix::trait_matrix_type& TraitMatrix::get_trait_matrix() const
{
    static trait_matrix_type traitmat(matrix.size());
    for (size_t i = 0; i < matrix.size(); ++i) {
        const trait_vector_type& vals = matrix[i].get_vals();
        traitmat[i].resize(vals.size());
        for (size_t j = 0; j < vals.size(); ++j) {
            traitmat[i][j] = vals[j];
        }
    }
    return(traitmat);
}


void TraitMatrix::realign_trait_matrix_to_leaf_ids() {
    assert(taxa != NULL);
    assert(taxa->taxa.size() > 0);
    const TaxonMatrix::taxon_id_vector_type& oldtid = get_taxon_ids();
    TaxonMatrix::taxon_id_vector_type newtid(oldtid.size(), TaxonMatrix::notset);
    for (size_t i = 0; i < newtid.size(); ++i) {
        newtid[i] = taxa->taxon_id_for_leaf_id(i);
        assert(newtid[i] != TaxonMatrix::notset);
    }
    realign_trait_matrix(newtid);
}


void TraitMatrix::realign_trait_matrix(
           const TaxonMatrix::taxon_id_vector_type& newtid)
{
    const TaxonMatrix::taxon_id_vector_type& oldtid = get_taxon_ids();
    assert(newtid.size() == oldtid.size());
    const int notset = -1;
    std::vector<int> newpos(oldtid.size(), notset);
    for (size_t opos = 0; opos < newpos.size(); ++opos) {
        for (size_t npos = 0; npos < newpos.size(); ++npos) {
            if (oldtid[opos] == newtid[npos]) {
                newpos[opos] = npos;
                break;
            }
        }
        if (newpos[opos] == notset) {
            std::cerr << "realign_trait_matrix: no newpos found for opos="
                      << opos << std::endl;
        }
        assert(newpos[opos] != notset);
    }

    trait_matrix_type mat(get_trait_matrix());
    transpose_trait_matrix(mat);
    trait_matrix_type newmat(mat.size());
    for (size_t i = 0; i < newpos.size(); ++i) {
        newmat[newpos[i]] = mat[i];
    }
    transpose_trait_matrix(newmat);
    set_trait_matrix(newmat);
    set_taxon_ids(newtid);
}


void TraitMatrix::set_trait_matrix(const TraitMatrix::trait_matrix_type& traitmat)
{
    assert(NULL != taxa);

    size_t num_entries = traitmat.size();
    size_t num_taxa = traitmat[0].size();
    if (taxa->taxa.size() != num_taxa) {
        std::cerr << "set_trait_matrix: no match: taxa->taxa.size()="
                  << taxa->taxa.size() << " and num_taxa=" << num_taxa
                  << std::endl;
        std::cerr << *taxa;
    }
    assert(taxa->taxa.size() == num_taxa);

    matrix.resize(num_entries);
    for (size_t i = 0; i < num_entries; ++i) {
        assert(traitmat[i].size() == num_taxa);
        matrix[i].set_id(i);
        matrix[i].vals.resize(num_taxa);
        matrix[i].init_vals(traitmat[i]);
    }
}


void TraitMatrix::transpose_trait_matrix(trait_matrix_type& traitmat)
{
    const size_t nrow = traitmat.size();
    const size_t ncol = traitmat[0].size();
    trait_matrix_type newmat(ncol);
    for (size_t i = 0; i < ncol; ++i) {
        newmat[i].resize(nrow);
        for (size_t j = 0; j < nrow; ++j) {
            newmat[i][j] = traitmat[j][i];
        }
    }
    traitmat = newmat;
}

#define DEBUG_READ_FILE 0

void TraitMatrix::write_file(const std::string& filename) const {
    std::ofstream file(filename.c_str());
    if (! file) {
        std::cerr << "TraitMatrix::write_file: could not open output file "
                  << filename << std::endl;
        assert(false);
    }
    
    const size_t num_taxa = taxon_ids.size();
    // line 1
    file << "traitsummary" << std::endl;
    // line 2
    for (size_t i = 0; i < num_taxa; ++i) {
        file << taxa->taxa[taxon_ids[i]].get_name();
        file << "\t";
    }
    file << "Count" << std::endl;
    // data lines
    for (size_t i = 0; i < matrix.size(); ++i) {
        const trait_vector_type& v = matrix[i].get_vals();
        for (size_t j = 0; j < num_taxa; ++j) {
            file << v[j] << "\t";
        }
        file << matrix[i].get_freq() << std::endl;
    }
}


void TraitMatrix::read_file(const std::string& filename)
{
    std::ifstream file(filename.c_str());
    if (! file) {
        std::cerr << "TraitMatrix::read_file: could not open input file "
                  << filename << std::endl;
        assert(false);
    }
    const int bufsz = 4096;
    char buf[bufsz+1];
    int numline = 0;
    int numtaxa = -1;
    int incident_count[2];
    incident_count[0] = 0;
    incident_count[1] = 0;
    const int maxfields = 40;
    bool traitsummary = false;
    bool headerseen = false;
    trait_matrix_type traitmat;
    std::vector< std::string > names;
    std::vector< freq_type > freqs;
    while (file.eof() == false) {
        ++numline;
        file.getline(buf, bufsz, '\n');
        std::istringstream inp(buf);

        std::vector<std::string> field(maxfields);
        int nf = 0;
        std::string w;
        while (inp >> w) {
            field[nf++] = w;
            if (nf >= maxfields)
                break;
        }
        if (nf == 0)
            continue;

        if (nf == 1 && field[0] == "traitsummary" && !traitsummary) {
            traitsummary = true;
            continue;
        }
        if (!traitsummary) {
            std::cerr << filename.c_str() << ":" << numline 
                      << " need initial `traitsummary' line" << std::endl;
            assert(false);
        }

        if (!headerseen) {
            // line contains taxa names, followed by "Count"
            numtaxa = nf - 1;
            for (int i = 0; i < numtaxa; ++i) {
                names.push_back(field[i]);
            }
            assert((int)names.size() == numtaxa);
            assert(field[numtaxa] == "Count" || 
                   field[numtaxa] == "Freqs" ||
                   field[numtaxa] == "Frequency");
            headerseen = true;
            continue;
        }

        trait_vector_type traitfield(numtaxa);
        freq_type freqfield;
        // peel off the fields one by one
        int local_incident_count[2];
        local_incident_count[0] = 0;
        local_incident_count[1] = 0;
        for (int i = 0; i < numtaxa; ++i) {
            std::istringstream inp(field[i].c_str());
            inp >> traitfield[i];
            ++local_incident_count[traitfield[i]];
        }
        traitmat.push_back(traitfield);

        std::istringstream inp2(field[numtaxa].c_str());
        inp2 >> freqfield;
        freqs.push_back(freqfield);
        incident_count[0] += local_incident_count[0]*freqfield;
        incident_count[1] += local_incident_count[1]*freqfield;
    }
    if (DEBUG_READ_FILE) {
        for (size_t i = 0; i < names.size(); ++i) {
            std::cout << names[i] << ",";
        }
        std::cout << std::endl;
        for (size_t i = 0; i < traitmat.size(); ++i) {
            std::cout << "traitmat[" << i << "]=";
            for (size_t j = 0; j < traitmat[i].size(); ++j) {
                std::cout << traitmat[i][j] << ",";
            }
            std::cout << " freq=" << freqs[i];
            std::cout << std::endl;
        }
    }

    fill(traitmat, freqs, names);
    realign_trait_matrix_to_leaf_ids();
    
    double sum = incident_count[0] + incident_count[1];
    std::cout << "TraitMatrix::read_file: Incident frequencies: sum="
              << long(sum) << " 0=" << incident_count[0]/sum
              << "  1=" << incident_count[1]/sum << std::endl;
}


void TraitMatrix::fill(const TraitMatrix::trait_matrix_type& traitmat,
                              const std::vector<freq_type>& freqs)
{
    assert(NULL != taxa);

    size_t num_entries = traitmat.size();
    assert(freqs.size() == num_entries);
    matrix.resize(0);  // to clear it, if it's not already
    matrix.resize(num_entries);
    
    for (size_t i = 0; i < num_entries; ++i) {
        matrix[i].set_freq(freqs[i]);
    }

    set_trait_matrix(traitmat);
    // realign_trait_matrix_to_leaf_ids();
}


void TraitMatrix::fill(const TraitMatrix::trait_matrix_type& traitmat,
                              const std::vector<freq_type>& freqs,
                              const std::vector<std::string>& names)
{
    // std::cout << "taxa=" << ((void*)taxa) << std::endl;
    // std::cout << taxa;
    assert(NULL != taxa);
    // assert(0 == taxa->taxa.size()); // may be able to soften this

    size_t num_taxa = traitmat[0].size();
    assert(names.size() == num_taxa);

    for (size_t i = 0; i < num_taxa; ++i) {
        TaxonMatrix::taxon_id_type tid = taxa->get_taxon_id(names[i]);
        if (tid < 0) {
            tid = taxa->add_taxon(names[i]);
        }
        taxon_ids[i] = tid;
    }

    fill(traitmat, freqs);
}


void TraitMatrix::set_taxa(TaxonMatrix& tm)
{
    taxa = &tm;
    for (size_t i = 0; i < taxa->taxa.size(); ++i) {
        if (taxa->taxa[i].get_leaf_id() == TaxonMatrix::notset) {
            std::cerr << "TraitMatrix::set_taxa(): taxa must have leaf_ids"
                      << " set, prior to this call, by calling "
                      << " PhyloTreeNode.find_taxa_for_leaves()" << std::endl;
        }
        assert(taxa->taxa[i].get_leaf_id() != TaxonMatrix::notset);
    }
    TaxonMatrix::taxon_id_vector_type v_tid;
    v_tid.assign(taxa->taxa.size(), TaxonMatrix::notset);
    for (size_t l = 0; l < v_tid.size(); ++l) {
        // loop through the columns from 0 -> taxon_ids.size()-1, looking
        // for a leaf_id that matches each one
        TaxonMatrix::taxon_id_type tid = taxa->taxon_id_for_leaf_id(l);
        if (tid == TaxonMatrix::notset) {
            std::cerr << "TraitMatris::set_taxa(): invalid taxon_id returned"
                      << " for taxa->taxon_id_for_leaf_id(" << l << ")"
                      << ", which may indicate taxa that are not also leaves"
                      << std::endl;
        }
        assert(tid != TaxonMatrix::notset);
        v_tid[l] = tid;
    }
    if (taxon_ids.size() > 0) {
        realign_trait_matrix(v_tid);
    } else {
        set_taxon_ids(v_tid);
    }
}


int TraitMatrix::map_taxon_id_to_col(int tid) const {
    assert(taxon_ids.size() > 0);
    for (size_t i = 0; i < taxon_ids.size(); ++i) {
        if (taxon_ids[i] == tid) {
            return((int)i);
        }
    }
    return(-1);  // out of bounds
}


void TraitMatrix::print(std::ostream& os) const {
    os << "max_num_states=" << max_num_states;
    os << " taxa=" << ((void*)taxa);
    os << " taxon_ids[0.." << (taxon_ids.size()-1) << "]={";
    for (size_t i = 0; i < taxon_ids.size(); ++i) {
        os << " " << taxon_ids[i];
    }
    os << "}";
    os << " leaf_ids[0.." << (taxon_ids.size()-1) << "]={";
    for (size_t i = 0; i < taxon_ids.size(); ++i) {
        os << " " << taxa->taxa[taxon_ids[i]].get_leaf_id();
    }
    os << "}";
    os << " valid_states[]={";
    for (int i = 0; i < max_num_states; ++i) {
        os << " " << valid_states[i];
    }
    os << "}";
    os << std::endl;
    for (size_t i = 0; i < matrix.size(); ++i) {
        os << "[" << i << "]";
        matrix[i].print(os);
    }
    os << std::endl;
}


std::ostream& operator<<(std::ostream& os, const TraitMatrix& tm) {
    os << "TraitMatrix:" << std::endl;
    tm.print(os);
    return(os);
}


void TraitMatrix::TraitMatrixEntry::print(std::ostream& os) const {
    os << "    id=" << id;
    for (size_t i = 0; i < vals.size(); ++i) {
        os << "   " << vals[i];
    }
    os << " freq=" << freq;
    os << std::endl;
}

