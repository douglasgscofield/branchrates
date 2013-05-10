#include "TaxonMatrix.h"

const int TaxonMatrix::notset = -1;

void TaxonMatrix::read_file(const std::string& filename)
{
    std::ifstream file(filename.c_str());
    assert(file);
    const int bufsz = 4096;
    char buf[bufsz+1];
    int numline = 0;
    const int maxfields = 2;
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
        switch(nf) {
            case 0: continue; break;
            case 2: add_taxon(field[0], field[1]); break;
            case 1: add_taxon(field[0]); break;
            default:
                std::cerr << filename.c_str() << ":" << numline 
                          << " requires 1 or 2 fields" << std::endl;
                assert(false);
                break;
        }
    }
}


std::ostream& operator<<(std::ostream& os, const TaxonMatrix& tm)
{
    os << "TaxonMatrix:" << std::endl;
    tm.print(os);
    return(os);
}


void TaxonMatrix::print(std::ostream& os) const
{
    for (int i = 0; i < taxa.size(); ++i) {
        taxa[i].print(os);
    }
    os << std::endl;
}

void TaxonMatrix::Taxon::print(std::ostream& os) const {
    os << " id=" << id;
    os << " name=:" << name << ":";
    os << " description=:" << description << ":";
    os << " leaf_id=" << leaf_id;
    os << std::endl;
}


