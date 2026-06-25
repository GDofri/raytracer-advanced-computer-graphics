//
// Created by dofri on 28.3.2024.
//


//
// Created by dofri on 15.3.2024.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/emitter.h>
#include <math.h>
#include <nori/bsdf.h>
#include <nori/sampler.h>
NORI_NAMESPACE_BEGIN

    /** MatsIntegrator, ...
     */
    class MatsIntegrator : public SampleIntegrator {
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {

            int grace = 3;
            Ray3f r = ray;
            Color3f totalRadiance = Color3f(0.0f, 0.0f, 0.0f);
            Color3f beta = Color3f(1.0f, 1.0f, 1.0f);

            for( int i = 0;;i++ )
            {
                /* Find the surface that is visible in the requested direction */

                Intersection its;
                if (!scene->rayIntersect(r, its))
                    return totalRadiance;

                if( its.mesh->isEmitter() )
                {
                    // We sampled a point on an area light so we add its emittance
                    if (its.toLocal(-r.d).z() > 0 )
                    {
                        totalRadiance += its.mesh->getEmitter()->getRadiance()*beta;
                    }
                }
                float contProb = 0.95;

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

                BSDFQueryRecord bRec(its.toLocal(-r.d));

                Color3f f = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

                Vector3f wo = bRec.wo;
                beta *= f;
                r = Ray3f(its.p, its.toWorld(wo), Epsilon, std::numeric_limits<float>::infinity());
            }
        }


    public:

        MatsIntegrator(const PropertyList &props) {
        }

        std::string toString() const override {
            return "MatsIntegrator[]";
        }


        std::string getName() const {
            return "MATS";
        }
    };

    NORI_REGISTER_CLASS(MatsIntegrator, "path_mats")
NORI_NAMESPACE_END