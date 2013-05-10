#include "BranchRateManager.h"


void BranchRateManager::set_unidir_rate(double r) {
    for (int i = 0; i < branchrates.size(); ++i) {
        set_unidir_rate(i, r);
    }
}


void BranchRateManager::set_unidir_rate(int brid, double r) {
    ratepv[branchrates[brid].get_forward_id()].set_ratep(r);
    ratepv[branchrates[brid].get_backward_id()].set_ratep(r);
}


void BranchRateManager::set_bidir_rates(double r_forw, double r_back) {
    for (int i = 0; i < branchrates.size(); ++i) {
        set_bidir_rates(i, r_forw, r_back);
    }
}


void BranchRateManager::set_bidir_rates(int brid, double r_forw,
                                               double r_back) {
    int fid = branchrates[brid].get_forward_id();
    int bid = branchrates[brid].get_backward_id();
    if (fid == bid) {
        std::cerr << "BranchRateManager::set_bidir_rate(): for branchrate_id="
                  << branchrates[brid].get_this_branchrate_id()
                  << ", only one ratep allocated" << std::endl;
    }
    assert(fid != bid);
    ratepv[fid].set_ratep(r_forw);
    ratepv[bid].set_ratep(r_back);
}


void BranchRateManager::calculate_p_transition() {
    for (int i = 0; i < branchrates.size(); ++i) {
        if (get_ratep_allocated(i) != true) {
            std::cerr << "not allocated: ratep for i=" << i << std::endl;
            print(std::cerr);
        }
        assert(get_ratep_allocated(i) == true);
        double f_rate = get_rate_forward(i);
        double b_rate = get_rate_backward(i);
        double rate_sum = f_rate + b_rate;
        double t1 = f_rate/rate_sum;
        double t2 = b_rate/rate_sum;
        double t22 = rate_sum * get_branchlength(i);
        double t3 = exp(-t22);
        double t4 = 1.0 - t3;
        double t5 = t1 * t4;
        double t6 = t2 * t4;
        double p01 = t5;
        double p00 = 1.0 - t5;
        double p10 = t6;
        double p11 = 1.0 - t6;
        set_p_transition(i, 0, 1, p01);
        set_p_transition(i, 0, 0, p00);
        set_p_transition(i, 1, 0, p10);
        set_p_transition(i, 1, 1, p11);
    }
}


void BranchRateManager::allocate_ratep_from_map(RatePMap& map)
{
    // a name-to-name map of relationships between branches and rate parameters
    if (ratepv.size() > 0) {
        std::cerr << "allocate_ratep_from_map: ratepv not empty" << std::endl;
    }
    assert(ratepv.size() == 0);
    for (int i = 0; i < map.size(); ++i) {
        //std::cerr << "map name " << i << " " << map.get_branchname(i) << std::endl;
        int brid;  // this will hold the branch referenced in map[i]
        for (brid = 0; brid < branchrates.size(); ++brid) {
            //std::cerr << "branch name " << brid << " " << branchrates[brid].get_name() << std::endl;
            if (branchrates[brid].get_name() == map.get_branchname(i)) {
                break;
            }
        }
        if (brid >= branchrates.size()) {
            std::cerr << "allocate_ratep_from_map: branch not found" << std::endl;
        }
        assert(brid < branchrates.size());

        // allocate the named rate parameter, or find it if already allocated
        int rid = ratepv.add_if_new_ratep(map.get_ratename(i));
        map.set_ratep_id(i, rid);

        if (map.get_is_forward(i) == true) {
            set_forward_id(brid, rid);
        } else if (map.get_is_backward(i) == true) {
            set_backward_id(brid, rid);
        } else {
            std::cerr << "allocate_ratep_from_map: unknown rate type" << std::endl;
        }
    }
    // go through all branchrates, make sure that forward and backward assigned
    for (int i = 0; i < branchrates.size(); ++i) {
        if (branchrates[i].get_forward_id() >= 0 &&
            branchrates[i].get_backward_id() >= 0) {

            branchrates[i].set_ratep_allocated(true);

        }
    }
}


void BranchRateManager::read_ratep_init_file(const std::string& filename)
{
    if (ratepv.size() == 0) {
        std::cerr << "BranchRateManager::read_ratep_init_file: rate parameters not yet allocated"
                  << std::endl;
        assert(false);
    }
    std::ifstream file(filename.c_str());
    if (! file) {
        std::cerr << "BranchRateManager::read_ratep_init_file: could not open input file "
                  << filename << std::endl;
        assert(false);
    }
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
        if (nf == 0)
            continue;
        if (nf != maxfields) {
            std::cerr << filename.c_str() << ":" << numline << ":"
                        << " need " << maxfields << " fields" << std::endl;
            assert(false);
        }
        double init_value;
        std::istringstream inp2(field[1].c_str());
        inp2 >> init_value;
        // so now field[0] is the rate parameter name, and init_value is
        // the value with which to initialize it
        int rid = ratepv.get_ratep_id(field[0]);
        if (rid < 0) {
            std::cerr << "unknown rate parameter name used for initialization: "
                      << field[0] << std::endl;
        }
        assert(rid >= 0);
        ratepv[rid].set_ratep(init_value);

        if (DEBUG_READ_RATEP_INIT_FILE) {
            // std::cout << "read_ratep_init_file: fields: ";
            // for (int i = 0; i < field.size(); ++i) {
            //     std::cout << field[i] << " ";
            // }
            // std::cout << std::endl;
            std::cout << "read_ratep_init_file: initialize " << field[0]
                      << " to " << init_value << std::endl;
        }
    }
    if (DEBUG_READ_RATEP_INIT_FILE) {
        std::cout << "Dump of RatePVector following initialization:" << std::endl;
        std::cout << ratepv;
        std::cout << std::endl;
    }
    
    ratepv.set_initialized(true);
}


void BranchRateManager::print(std::ostream& os) const {
    os << "name_prefix=:" << name_prefix << ":";
    os << " name_forward=:" << name_forward << ":";
    os << " name_backward=:" << name_backward << ":";
    os << std::endl;
    for (int i = 0; i < branchrates.size(); ++i) {
        print(branchrates[i].get_this_branchrate_id(), os);
        os << std::endl;
    }
}


void BranchRateManager::print(int brid, std::ostream& os) const {
    static const int annotate_ratep_entries = 1;
    if (brid < 0) {
        os << "<branchrate_id not set>";
        return;
    }
    branchrates[brid].print(os);
    if (annotate_ratep_entries &&
        branchrates[brid].get_ratep_allocated() == true) {
        int id = branchrates[brid].get_forward_id();
        os << " RatePEntry[" << id << "]={";
        //((const RatePVector::RatePEntry&)ratepv[id]).print(os);
        (((BranchRateManager*)this)->ratepv)[id].print(os);
        os << "}";
        id = branchrates[brid].get_backward_id();
        os << " RatePEntry[" << id << "]={";
        (((BranchRateManager*)this)->ratepv)[id].print(os);
        //((const RatePVector::RatePEntry&)ratepv[id]).print(os);
        os << "}";
    }
}


std::ostream& operator<<(std::ostream& os, const BranchRateManager& brm)
{
    os << "BranchRateManager: ";
    brm.print(os);
    os << std::endl;
}


