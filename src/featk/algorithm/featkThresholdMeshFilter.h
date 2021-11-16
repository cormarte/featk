#ifndef FEATKTHRESHOLDMESHFILTER_H
#define FEATKTHRESHOLDMESHFILTER_H

#include <featk/algorithm/featkMeshConsumerBase.h>
#include <featk/algorithm/featkMeshProducerBase.h>
#include <featk/core/featkGlobal.h>
#include <featk/geometry/featkMesh.h>

template<unsigned int Dimension>
class featkThresholdMeshFilter : public virtual featkMeshConsumerBase<Dimension>, public virtual featkMeshProducerBase<Dimension> {

    public:

        featkThresholdMeshFilter();
        ~featkThresholdMeshFilter();

        FEATK_EXPORT void execute();

        std::vector<std::map<size_t, size_t>> getElementMaps();
        std::vector<std::map<size_t, size_t>> getNodeMaps();
        std::vector<std::map<size_t, size_t>> getReversedElementMaps();
        std::vector<std::map<size_t, size_t>> getReversedNodeMaps();
        void setAttributeName(std::string name);
        void setThreshold(double threshold);

    private:

        std::string attributeName;
        std::vector<std::map<size_t, size_t>> elementMaps;
        std::vector<std::map<size_t, size_t>> nodeMaps;
        std::vector<std::map<size_t, size_t>> reversedElementMaps;
        std::vector<std::map<size_t, size_t>> reversedNodeMaps;
        double threshold;
};

template<unsigned int Dimension>
featkThresholdMeshFilter<Dimension>::featkThresholdMeshFilter() {

    this->attributeName = "";
    this->threshold = 0.0;

    this->setNumberOfOutputMeshes(2);
}

template<unsigned int Dimension>
featkThresholdMeshFilter<Dimension>::~featkThresholdMeshFilter() {

}

template<unsigned int Dimension>
std::vector<std::map<size_t, size_t>> featkThresholdMeshFilter<Dimension>::getElementMaps() {

    return this->elementMaps;
}

template<unsigned int Dimension>
std::vector<std::map<size_t, size_t>> featkThresholdMeshFilter<Dimension>::getNodeMaps() {

    return this->nodeMaps;
}

template<unsigned int Dimension>
std::vector<std::map<size_t, size_t>> featkThresholdMeshFilter<Dimension>::getReversedElementMaps() {

    return this->reversedElementMaps;
}

template<unsigned int Dimension>
std::vector<std::map<size_t, size_t>> featkThresholdMeshFilter<Dimension>::getReversedNodeMaps() {

    return this->reversedNodeMaps;
}

template<unsigned int Dimension>
void featkThresholdMeshFilter<Dimension>::setAttributeName(std::string name) {

    this->attributeName = name;
}

template<unsigned int Dimension>
void featkThresholdMeshFilter<Dimension>::setThreshold(double threshold) {

    this->threshold = threshold;
}

#endif // FEATKTHRESHOLDMESHFILTER_H
