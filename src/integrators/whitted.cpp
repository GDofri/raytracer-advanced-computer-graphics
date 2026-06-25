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

    /** WhittedIntegrator, ...
     */
    class WhittedIntegrator : public SampleIntegrator {
    public:
        WhittedIntegrator
        (const PropertyList &props) {
        }

        // Todo It would probably be better to do this as a loop as opposed to the current recursive implementation.
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const override {
            /* Find the surface that is visible in the requested direction */
            Color3f totalRadiance = Color3f(0.0f, 0.0f, 0.0f);
            Intersection its;
            if (!scene->rayIntersect(ray, its))
                return totalRadiance;

            auto [sampledEmitter, emitterPdf ] = scene->sampleEmitter(sampler->next1D());

            if( its.mesh->isEmitter() )
            {
                // We sampled a point on an area light so we add its emittance
                if (its.toLocal(-ray.d).z() > 0 )
                {
                    totalRadiance += its.mesh->getEmitter()->getRadiance();
                }
            }

            if( its.mesh->getBSDF()->isDiffuse() )
            {
                // Check if light source is occluded from the intersection point.
                EmitterQueryRecord eqr;
                sampledEmitter->getEmitter()->sample( eqr, *sampler);
                float lightLocPdf = eqr.pdfLoc;
                Vector3f vecToLight = eqr.loc - its.p;
                Vector3f dirToLight = vecToLight.normalized();
                Ray3f shadowRay = Ray3f(its.p, dirToLight, Epsilon, vecToLight.norm() - Epsilon);


                if (scene->rayIntersect(shadowRay))
                {
                    // Light is occluded
                    return totalRadiance;
                }


                // Light is not occluded
                float cosI = its.shFrame.n.dot(dirToLight);
                float cosO = eqr.normal.dot(-dirToLight);
                float distToLightSquared = (eqr.loc - its.p).squaredNorm();

                float geometricTerm = fmax(0.0f, cosI * cosO) / distToLightSquared;

                BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(dirToLight), ESolidAngle);

                totalRadiance += its.mesh->getBSDF()->eval(bRec) *
                                 geometricTerm *
                                 eqr.radiance /
                        ( emitterPdf * lightLocPdf );

                return totalRadiance;
            }

            // BSDF is specular.
            BSDFQueryRecord bRec(its.toLocal(-ray.d));
            Color3f bsdfColor = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            float continueProb = 0.95;
            if( sampler->next1D() < continueProb )
            {
                return (1.0f/(continueProb))*Li(scene, sampler, Ray3f(its.p, its.toWorld(bRec.wo)));
            }
            else
            {
                return 0;
            }
        }

        std::string toString() const override {
            return "WhittedIntegrator[]";
        }


        std::string getName() const {
            return "Whitted";
        }
    };

    NORI_REGISTER_CLASS(WhittedIntegrator, "whitted")
NORI_NAMESPACE_END