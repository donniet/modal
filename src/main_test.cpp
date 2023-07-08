#include "modal.hpp"

#include <gtest/gtest.h>

#include <utility>
#include <cmath>
#include <algorithm>
#include <execution>

using namespace modal;

TEST(guassian, guass) {
    std::cout << "invsqrt2" << invsqrt2pi << std::endl;

    auto a = gauss(0., 1.);

    EXPECT_EQ(a, 1. / std::sqrt(2. * pi));
}

TEST(gaussian, instancing) {
    normal_distribution<double> a(0.,1.), b(0.,1.), c(1., 1.);

    auto y = a.value(0);

    EXPECT_EQ(y, invsqrt2pi);

    normal_distribution<double> d;
    double error;
    std::tie(d, error) = mix(a, b);

    EXPECT_EQ(error, 0.);
}

TEST(gaussian, inserting) {
    modal::modal m;
    
    std::random_device rd;
    std::normal_distribution<> d{1., 1.};

    double mean = 0.;
    double m2 = 0.;

    const int count = 15;

    for(int i = 0; i < count; i++) {
        double x = d(rd);
        double delta = x - mean;
        mean += delta / (double)(i+1);
        m2 += delta*(x-mean);

        std::cout << mean << " " << m2 << std::endl;

        m.insert(x);
    }

    std::cout << "actual mean: " << mean << std::endl;
    std::cout << "stddev:      " << sqrt(m2 / (count-1)) << std::endl;

    std::cout << "mean:        " << m.root->mean() << std::endl;
    std::cout << "stddev:      " << m.root->standard_deviation() << std::endl;
    std::cout << "error:       " << m.root->mixture_error() << std::endl;

    // m.visit([](int depth, double mean, double stdev, double scale, double err)->bool {
    //     for(; depth > 0; depth--) std::cout << " ";
    //     std::cout << "[" << mean << " ; " << stdev << " x " << scale << " ^ " << err << "]" << std::endl;

    //     return true;
    // });

}

TEST(gaussian, extracting) {
    modal::modal m;
    
    std::random_device rd;
    std::normal_distribution<> d{2., 1.};
    std::normal_distribution<> e{-2., 1.};
    std::uniform_int_distribution<> c{0, 1};

    const int count = 150000;

    std::vector<double> data(count);

    // std::generate_n(std::execution::par_unseq, data.begin(), count, [&c, &d, &e, &rd]() -> double {
    //     if(c(rd) == 0) {
    //         return d(rd);
    //     }
    //     return e(rd);
    // });

    // std::for_each(data.begin(), data.end(), [&m](double x) {
    //     m.insert(x);
    // });

    for(int i = 0; i < count; i++) {
        double x;
        if(c(rd) == 0) {
            x = d(rd);
        } else {
            x = e(rd);
        }

        m.insert(x);
    }

    auto modes = m.extract();

    for(auto const & d : modes) {
        std::cout << "mean: " << d.mean << " sdev: " << d.standard_deviation << " scale: " << d.scale << std::endl;
    }

    int parentLess = 0, parentMore = 0, mixed = 0, leaf = 0;
    m.visit_nodes([&](int depth, node * here) -> bool {
        if(here->l == nullptr) {
            leaf++;
        } else if(here->mixture_error() < here->l->mixture_error()) {
            if(here->mixture_error() < here->r->mixture_error()) {
                parentLess++;
            } else {
                mixed++;
            }
        } else if(here->mixture_error() < here->r->mixture_error()) {
            mixed++;
        } else {
            parentMore++;
        }
        return true;
    });

    std::cout << "parentLess: " << parentLess << std::endl
              << "parentMore: " << parentMore << std::endl
              << "mixed:      " << mixed << std::endl
              << "leaf:       " << leaf << std::endl;

    // m.visit([](int depth, double mean, double stdev, double scale, double err) -> bool {
    //     for(; depth > 0; depth--) std::cout << " ";
    //     std::cout << "[" << mean << " ; " << stdev << " x " << scale << " ^ " << err << "]" << std::endl;

    //     return depth < 7;
    // });

}