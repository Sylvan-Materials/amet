#pragma once
#include <functional>

namespace amet{

    template<typename G>
    class vf2 {
        public:
        vf2() = default;
        vf2(const vf2& orig) = delete;
        virtual ~fv2() = default;

        void operator()(std::function<void()> apply){

        }

    };
}