#ifndef RATEPVECTOR_H
#define RATEPVECTOR_H

#include <vector>
#include <iostream>
#include <string>
#include <cassert>

class RatePVector {
    private:
        std::string                   name_prefix;

    public:
        typedef std::vector<double>   v_ratep_type;
        typedef std::vector<v_ratep_type>
                                      m_ratep_type;
        class RatePEntry {
            private:
                int                   ratep_id;
                double                ratep;
                double                prev_ratep;
                std::string           name;
                std::string           description;

            public:
                int                   get_ratep_id() const;
                void                  set_ratep_id(const int rid);
                double                get_ratep() const;
                void                  set_ratep(const double arg);
                void                  update_ratep(const double delta);
                double                get_prev_ratep() const;
                void                  set_prev_ratep(const double arg);
                const std::string&    get_name() const;
                void                  set_name(const std::string& n);
                const std::string&    get_description() const;
                void                  set_description(const std::string& d);
                void                  print(std::ostream& = std::cout) const;
        };
        typedef std::vector<RatePEntry>
                                      ratepvector_type;

    private:
        ratepvector_type              ratepvec;
        bool                          prev_in_use;
        bool                          initialized;

    protected:
        void                          store();
        void                          restore();

    public:
        RatePVector();

        int                      get_ratep_id(const std::string& n) const;
        void                     set_initialized(bool arg = true);
        bool                     get_initialized();
        int                      add_new_ratep(std::string n = "",
                                               std::string d = "",
                                               double val = 0.0);
        int                      add_if_new_ratep(std::string n = "",
                                                  std::string d = "",
                                                  double val = 0.0);
        int                      size() const;
        void                     assign(double arg);
        void                     assign(const v_ratep_type& v);
        void                     update(const v_ratep_type& v);
        const v_ratep_type&      get_v_ratep() const;

        void                     print(std::ostream& os = std::cout) const;
        RatePEntry&              operator[](int index);
        RatePVector&             operator=(const RatePVector& rpv);
        friend std::ostream&     operator<<(std::ostream& os,
                                            const RatePVector& v);

};

////////////////////////////
//////////////////////////// RatePVector::RatePEntry methods
////////////////////////////

inline void RatePVector::RatePEntry::print(std::ostream& os) const {
    os << "ratep_id=" << ratep_id;
    os << " ratep=" << ratep;
    os << " prev_ratep=" << prev_ratep;
    os << " name=:" << name << ":";
    os << " description=:" << description << ":";
}

inline int RatePVector::RatePEntry::get_ratep_id() const {
    return(ratep_id);
}

inline void RatePVector::RatePEntry::set_ratep_id(const int rid) {
    ratep_id = rid;
}

inline double RatePVector::RatePEntry::get_ratep() const {
    return(ratep);
}

inline void RatePVector::RatePEntry::set_ratep(const double arg) {
    ratep = arg;
}

inline void RatePVector::RatePEntry::update_ratep(const double delta) {
    set_ratep(get_ratep() + delta);
}

inline double RatePVector::RatePEntry::get_prev_ratep() const {
    return(prev_ratep);
}

inline void RatePVector::RatePEntry::set_prev_ratep(const double arg) {
    prev_ratep = arg;
}

inline const std::string& RatePVector::RatePEntry::get_name() const {
    return(name);
}

inline void RatePVector::RatePEntry::set_name(const std::string& n) {
    name = n;
}

inline const std::string& RatePVector::RatePEntry::get_description() const {
    return(description);
}

inline void RatePVector::RatePEntry::set_description(const std::string& d) {
    description = d;
}

////////////////////////////
//////////////////////////// RatePVector methods
////////////////////////////

inline RatePVector::RatePEntry& RatePVector::operator[](int index) {
    return(ratepvec[index]);
}

inline RatePVector& RatePVector::operator=(const RatePVector& rpv) {
    if (&rpv != this) {
        ratepvec = rpv.ratepvec;
        prev_in_use = rpv.prev_in_use;
    }
    return(*this);
}

inline void RatePVector::set_initialized(bool arg) {
    initialized = arg;
}

inline bool RatePVector::get_initialized() {
    return(initialized);
}

inline int RatePVector::get_ratep_id(const std::string& n) const {
    for (int i = 0; i < ratepvec.size(); ++i) {
        if (n == ratepvec[i].get_name()) {
            return(ratepvec[i].get_ratep_id());
        }
    }
    return(-1);
}

inline int RatePVector::add_if_new_ratep(std::string n, std::string d,
                                         double val) {
    int ratep_id = get_ratep_id(n);
    if (ratep_id >= 0) {
        return(ratep_id);
    }
    return(add_new_ratep(n, d, val));
}

inline int RatePVector::add_new_ratep(std::string n, std::string d, double val) {
    if (get_ratep_id(n) >= 0) {
        std::cerr << "add_new_ratep(): duplicate name :" << n
                  << ":" << std::endl;
    }
    assert(get_ratep_id(n) < 0);  // not a duplicate entry
    int end = ratepvec.size();
    ratepvec.resize(end + 1);
    ratepvec[end].set_ratep_id(end);
    ratepvec[end].set_name(n);
    ratepvec[end].set_description(d);
    ratepvec[end].set_ratep(val);
    ratepvec[end].set_prev_ratep(0.0);
    return(end);
}

inline int RatePVector::size() const {
    return(ratepvec.size());
}

inline void RatePVector::assign(double arg) {
    for (int i = 0; i < ratepvec.size(); ++i) {
        ratepvec[i].set_ratep(arg);
    }
}

inline const RatePVector::v_ratep_type& RatePVector::get_v_ratep() const {
    static v_ratep_type v;
    v.resize(size());
    for (int i = 0; i < size(); ++i) {
        v[i] = ratepvec[i].get_ratep();
    }
    return(v);
}

inline void RatePVector::assign(const v_ratep_type& v) {
    assert(v.size() == size());
    for (int i = 0; i < size(); ++i) {
        ratepvec[i].set_ratep(v[i]);
    }
}

inline void RatePVector::update(const v_ratep_type& v) {
    assert(v.size() == size());
    for (int i = 0; i < size(); ++i) {
        ratepvec[i].update_ratep(v[i]);
    }
}

inline RatePVector::RatePVector() : name_prefix("rate_"), prev_in_use(false),
    initialized(false)
{
    /* EMPTY */
}

inline void RatePVector::store() {
    for (int i = 0; i < size(); ++i) {
        ratepvec[i].set_prev_ratep(ratepvec[i].get_ratep());
    }
    prev_in_use = true;
}

inline void RatePVector::restore() {
    assert(prev_in_use == true);
    for (int i = 0; i < ratepvec.size(); ++i) {
        ratepvec[i].set_ratep(ratepvec[i].get_prev_ratep());
    }
    prev_in_use = false;
}

#endif // RATEVECTOR_H

