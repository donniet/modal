#pragma once

#include <cmath>
#include <numbers>
#include <utility>
#include <tuple>
#include <random>
#include <vector>
#include <limits>

namespace modal {

using std::numbers::pi;
using std::numbers::sqrt2;
using std::numbers::e;

using std::sqrt;
using std::exp;
using std::erf;
using std::pair;
using std::tie;
using std::vector;
using std::make_pair;

constexpr double sqrtpi = sqrt(pi);
constexpr double sqrt2pi = sqrt(2. * pi);
constexpr double invsqrt2pi = 1. / sqrt2pi;

template<typename T>
inline T sq(T x) {
    return x * x;
}

template<typename T>
inline T gauss(T x) {
    return invsqrt2pi * exp(-0.5 * x * x);
}
template<typename T>
inline T gauss(T x, T s) {
    x /= s;
    return invsqrt2pi * exp(-0.5 * x * x) / s;
}

template<typename T>
class normal_distribution {
public:
    T mean;
    T standard_deviation;
    T scale;
public:
    normal_distribution() 
        : mean(0), standard_deviation(0), scale(1.) 
    { }

    normal_distribution(T mean, T standard_deviation)
        : mean(mean), standard_deviation(standard_deviation), scale(1.) 
    { }

    normal_distribution(T mean, T standard_deviation, T scale)
        : mean(mean), standard_deviation(standard_deviation), scale(scale) 
    { }

    inline T value(T x) const {
        if(standard_deviation == 0.) {
            if(x == mean) return INFINITY;
            return 0;
        }
        return scale * gauss(x-mean, standard_deviation); 
    }

    inline T probability(T x) const {
        if(standard_deviation == 0.) return 0.;
        return std::erf((x - mean) / sqrt2 / standard_deviation);
    }

    inline T probability1(T x, double alpha=1.) {
        if(standard_deviation == 0.) return 1.;
        
        // adjust the standard deviation by the scale
        standard_deviation = alpha * scale / (alpha * scale - 1.);

        return std::erf((x - mean) / sqrt2 / standard_deviation);

    }
};

template<typename T>
pair<normal_distribution<T>, T> 
mix(normal_distribution<T> const & l, normal_distribution<T> const & r) 
{
    T tot = l.scale + r.scale;
    T m = (l.scale * l.mean + r.scale * r.mean) / tot;
    T v =   (l.scale * (sq(l.standard_deviation) + sq(m - l.mean)) + 
             r.scale * (sq(r.standard_deviation) + sq(m - r.mean)))/tot;
    T s = sqrt(v);
    auto ret = normal_distribution<T>(m, s, tot);

    // calculate the error in the mixture
    double error = 0.;


    if(l.standard_deviation == 0. || r.standard_deviation == 0.) {
        error = std::numeric_limits<T>::max();
    } else {
        error += sq(l.scale / tot) / 2. / sqrtpi / l.standard_deviation;
        error += sq(r.scale / tot) / 2. / sqrtpi / r.standard_deviation;
        error += 1. / 2. / sqrtpi / s;
        error += 2. * l.scale * r.scale / sq(tot) * gauss(l.mean - r.mean, sqrt(sq(l.standard_deviation) + sq(r.standard_deviation)));
        error -= 2. * l.scale / tot * gauss(m - l.mean, sqrt(v + sq(l.standard_deviation)));
        error -= 2. * r.scale / tot * gauss(m - r.mean, sqrt(v + sq(r.standard_deviation)));
    }

    return std::make_pair(ret, error);
}

double mixture_error(double m, double v, double alpha, double m1, double v1, double m2, double v2) {
    double error = 0.;
    if(v1 == 0. || v2 == 0.) return 0; //std::numeric_limits<double>::max();

    double s = sqrt(v);
    double s1 = sqrt(v1);
    double s2 = sqrt(v2);

    error =   sq(alpha)      / 2. / sqrtpi / s1
            + sq(1. - alpha) / 2. / sqrtpi / s2
            + 1.             / 2. / sqrtpi / s
            + 2. * alpha * (1. - alpha) * gauss(m1 - m2, sqrt(v1 + v2))
            - 2. * alpha                * gauss(m  - m1, sqrt(v  + v1))
            - 2.         * (1. - alpha) * gauss(m  - m2, sqrt(v  + v2));

    return error;
}

class modal;
class node;


class node {
    friend class modal;
public:
    double m;
    double m2;
    double count;
    double err;

    node * l, * r;
public:
    node(node * left, node * right)
        :  m(0.), m2(0.), count(0), l(left), r(right)
    {
        if(l == nullptr && r == nullptr) {
            return;
        }
        
        count = l->count + r->count;
        m = (l->count * l->m + r->count * r->m) / count;
        m2 = l->m2 + r->m2 + l->count * r->count * sq(l->m - r->m) / count;
        err = ::modal::mixture_error(m, variance(), l->count / count, l->m, l->variance(), r->m, r->variance());
    }

    node(double sample) 
        : m(sample), m2(0), count(1), l(nullptr), r(nullptr), err(std::numeric_limits<double>::max())
    { }

    double mean() const { return m; }
    double variance() const { 
        if(count < 2) return 0;

        return m2 / (count-1);
    }
    double standard_deviation() const { return sqrt(variance()); }
    double mixture_error() const { return err; }

    inline double probability1(double x, double alpha=1.) const {
        auto sdev = standard_deviation();
        if(sdev == 0.) return 1.;
        
        // adjust the standard deviation by the scale
        sdev = alpha * count / (alpha * count - 1.);

        return std::erf((x - mean()) / sqrt2 / sdev);
    }

    inline double density(double x) const {
        return gauss(x - mean(), standard_deviation());
    }
};


class modal {
public:
    node * root;

public:
    modal() : root(nullptr) { }
    ~modal() {
        delete_helper(&root);
    }

    void insert(double sample) {
        insert_helper(&root, sample);
       
    }

    template<typename Visitor>
    void visit(Visitor v) {
        visit_helper(0, root, v);
    }

    template<typename Visitor>
    void visit_nodes(Visitor v) {
        std::vector<pair<int, node*>> stack;
        stack.push_back(make_pair(0, root));
        int cur;
        node * top;

        while(stack.size() > 0) {
            tie(cur, top) = stack.back();
            stack.pop_back();
            
            if(!v(cur, top)) {
                break;
            }

            if(top->l != nullptr) stack.push_back(make_pair(cur+1, top->l));
            if(top->r != nullptr) stack.push_back(make_pair(cur+1, top->r));
        }
    }

    vector<normal_distribution<double>> extract() const {
        vector<normal_distribution<double>> ret;

        extract_helper(root, ret);
        
        return ret;
    }
private:
    void extract_helper(node * here, vector<normal_distribution<double>> & ret) const {
        if(here == nullptr) return;

        // we will extract if our error is less than both our child errors
        if(here->l == nullptr || here->r == nullptr || 
            (here->mixture_error() < here->l->mixture_error() && here->mixture_error() < here->r->mixture_error())) 
        {
            ret.push_back(normal_distribution<double>(here->mean(), here->standard_deviation(), here->count));
            return;
        }

        extract_helper(here->l, ret);
        extract_helper(here->r, ret);
    }

    template<typename Visitor>
    void visit_helper(int depth, node * here, Visitor v) {
        if(here == nullptr) return;

        bool cont = v(depth, here->mean(), here->standard_deviation(), here->count, here->mixture_error());

        if(!cont) return;

        visit_helper(depth+1, here->l, v);
        visit_helper(depth+1, here->r, v);
    }
    void delete_helper(node ** here) {
        if(*here == nullptr) {
            return;
        }

        delete_helper(&(*here)->l);
        delete_helper(&(*here)->r);

        delete *here;
        *here = nullptr;
    }
    void insert_helper(node ** here, double sample) {
        node * h = *here;

        if(h == nullptr) {
            *here = new node(sample);
            return;
        }

        if(h->l == nullptr && h->r == nullptr) {
            auto n = new node(sample);
            auto m = new node(h, n);

            *here = m;
            h = *here;
        } else {
            h->count += 1;
            double delta = sample - h->m;
            h->m += delta / h->count;
            h->m2 += delta*(sample-h->m);

            double left = h->l->count * h->l->probability1(sample);
            double right = h->r->count * h->r->probability1(sample);

            std::minstd_rand gen;
            std::uniform_real_distribution<> dis(0, left + right);

            if(dis(gen) < left) {
                insert_helper(&(h->l), sample);
            } else {
                insert_helper(&(h->r), sample);
            }

            h->err = ::modal::mixture_error(h->mean(), h->variance(), h->l->count / h->count, h->l->mean(), h->l->variance(), h->r->mean(), h->r->variance());

            // handle the mixed states
            if(h->err < h->l->err && h->err > h->r->err) {

            } else if(h->err > h->l->err && h->err < h->r->err) {

            }
        }
    }
};




} // modal