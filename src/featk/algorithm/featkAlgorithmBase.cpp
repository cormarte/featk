#include <featk/algorithm/featkAlgorithmBase.h>

#include <iostream>

featkAlgorithmBase::featkAlgorithmBase() {

}

featkAlgorithmBase::~featkAlgorithmBase() {

}

bool featkAlgorithmBase::check() {

    std::cout << "featkAlgorithmBase: Warning: check() method not reimplemented." << std::endl;

    return true;
}

bool featkAlgorithmBase::update() {

    bool ret = false;

    if (this->check()) {

        this->execute();
        ret = true;
    }

    return ret;
}
