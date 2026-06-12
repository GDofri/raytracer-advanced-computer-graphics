//
// Created by dofri on 28.3.2024.
//

#include "nori/area.h"
#include "nori/warp.h"
NORI_NAMESPACE_BEGIN

AreaLight::AreaLight(const PropertyList& propList)
{
    m_radiance = propList.getColor("radiance");
}

void AreaLight::sample(EmitterQueryRecord &eqr, Sampler &sampler)
{
    auto [loc, normal, pdf] = m_mesh->sample(sampler);
    eqr = {loc, normal, eval(eqr), Vector3f::Zero(), pdf};
}

Color3f AreaLight::getRadiance() const
{
    return m_radiance;
}

float AreaLight::pdfPos(const EmitterQueryRecord &eqr) const
{
    return m_mesh->pdf();
}

float AreaLight::pdfDir(const EmitterQueryRecord &eqr) const
{
    return Warp::squareToCosineHemispherePdf( Frame(eqr.normal).toLocal(eqr.woWorld));
}

Color3f AreaLight::eval(const EmitterQueryRecord &record) const
{
    return m_radiance;
}

void AreaLight::sampleRay(EmitterQueryRecord &eqr, Sampler &sample, Ray3f &ray)
{
    auto [loc, normal, pdfLoc] = m_mesh->sample(sample);

    Frame shFrame = Frame(normal);

    Vector3f dir = Warp::squareToCosineHemisphere(sample.next2D());
    float pdfDir =  Warp::squareToCosineHemispherePdf(dir);
    Vector3f dirWorld = shFrame.toWorld(dir);
    eqr = {loc, normal, eval(eqr), dirWorld, pdfLoc, pdfDir };

    ray = Ray3f(loc, dirWorld, Epsilon, std::numeric_limits<float>::infinity());
}


NORI_NAMESPACE_END