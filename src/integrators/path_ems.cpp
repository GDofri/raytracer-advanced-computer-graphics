//
// Created by dofri on 28.3.2024.
//


//
// Created by dofri on 15.3.2024.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
NORI_NAMESPACE_BEGIN

    /** EmsIntegrator, ...
     */
    class EMSIntegrator : public SampleIntegrator {
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {

            int grace = 3;
            auto r = Ray3f(ray);
            float eta = 1.0f;
            Color3f totalRadiance = Color3f(0.0f, 0.0f, 0.0f);
            Color3f beta = Color3f(1.0f, 1.0f, 1.0f);
            bool lastWasSpecular = true;
            int maxDepth = 2;
            bool useMaxDepth = false;

            for( int i = 0; i < maxDepth || !useMaxDepth; i++ )
            {
                /* Find the surface that is visible in the requested direction */

                Intersection its;
                if (!scene->rayIntersect(r, its))
                    return totalRadiance;

                if( lastWasSpecular && its.mesh->isEmitter() )
                {
                    if (its.toLocal(-r.d).z() > 0 )
                    {
                        totalRadiance += its.mesh->getEmitter()->getRadiance()*beta;
                    }
                }

                // Direct lighting
                if( its.mesh->getBSDF()->isDiffuse()/* && its.mesh->getEmitter() == nullptr*/)
                {
                    auto [emitter, emitterPdf] = scene->sampleEmitter(sampler->next1D());
                    EmitterQueryRecord eqr;
                    emitter->getEmitter()->sample( eqr, *sampler);
                    float lightLocPdf = eqr.pdfLoc;
                    Vector3f vecToLight = eqr.loc - its.p;
                    Vector3f dirToLight = vecToLight.normalized();
                    Ray3f shadowRay = Ray3f(its.p, dirToLight, Epsilon, vecToLight.norm() - Epsilon);
                    if (!scene->rayIntersect(shadowRay))
                    {
                        // Light is not occluded
                        float cosI = its.shFrame.n.dot(dirToLight);
                        float cosO = eqr.normal.dot(-dirToLight);
                        float distToLightSquared = (eqr.loc - its.p).squaredNorm();
                        float geometricTerm = fmax(0.0f, cosI * cosO) / distToLightSquared;
                        BSDFQueryRecord bRec(its.toLocal(-r.d), its.toLocal(dirToLight), ESolidAngle);
                        totalRadiance += its.mesh->getBSDF()->eval(bRec) *
                                         beta *
                                         geometricTerm *
                                         eqr.radiance /
                                       ( emitterPdf * lightLocPdf);
                    }
                    lastWasSpecular = false;
                }
                else
                {
                    lastWasSpecular = true;
                }

                // Russian roulette
                float contProb = std::min(beta.maxCoeff() * eta * eta ,0.99f);
                if(i >= grace)
                {
                    if(sampler->next1D() >= contProb)
                    {
                        return totalRadiance;
                    }
                    else
                    {
                        beta /= contProb;
                    }
                }

                // Sample the BSDF and update the path
                BSDFQueryRecord bRec(its.toLocal(-r.d));

                Color3f f = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                Vector3f wo = bRec.wo;
                eta /= bRec.eta;
                beta *= f;
                r = Ray3f(its.p, its.toWorld(wo));

            }
            return totalRadiance;
        }

    public:

        EMSIntegrator(const PropertyList &props) {
        }

        std::string toString() const override {
            return "EmsIntegrator[]";
        }

        std::string getName() const {
            return "EMS";
        }
    };

    NORI_REGISTER_CLASS(EMSIntegrator, "path_ems")
NORI_NAMESPACE_END
