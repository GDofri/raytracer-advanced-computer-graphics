//
// Created by dofri on 28.3.2024.
//

#ifndef NORI_AREA_H
#define NORI_AREA_H


#include "emitter.h"

NORI_NAMESPACE_BEGIN

/**
 * Area light
 */
class AreaLight : public Emitter
{

public:

    /**
     * Area light constructor
     * @param propList
     */
    AreaLight(const PropertyList& propList);

    float pdfPos(const EmitterQueryRecord &eqr) const override;

    float pdfDir(const EmitterQueryRecord &eqr) const override;

    /**
     * Sample the emitter..
     * @param sampler
     * @return EmitterQueryRecord with the location, normal, radiance and the pdf value of the sample.
     */
    void sample(EmitterQueryRecord &eqr, Sampler &sampler) override;

    /**
     * Sample an outgoing ray from the emitter
     * @param eqr
     * @param sampler
     * @param ray
     */
    void sampleRay(EmitterQueryRecord &eqr, Sampler &sampler, Ray3f &ray);


    Color3f eval(const EmitterQueryRecord &record) const override;

    /**
     * Get the radiance of the emitter.
     * @return The radiance of the emitter.
     */
    Color3f getRadiance() const override;

    std::string toString() const override
    {
        return tfm::format( "AreaLight[radiance=%s]", m_radiance.toString());
    }
private:

    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area")

NORI_NAMESPACE_END

#endif //NORI_AREA_H
