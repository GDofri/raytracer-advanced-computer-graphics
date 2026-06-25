
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

    /** MisIntigrator, ...
     */
    class MisIntigrator : public SampleIntegrator {
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &initialRay) const override {


            int grace = 5;
            Ray3f ray = Ray3f(initialRay);
            float eta = 1.0f;
            Color3f totalRadiance = Color3f(0.0f, 0.0f, 0.0f);
            Color3f beta = Color3f(1.0f, 1.0f, 1.0f);
            bool lastWasSpecular;

            // Declare MIS weight for brdf here because it needs
            // to persist between iterations.
            float wBRDF = 1.0f;

            Intersection its;
            if (!scene->rayIntersect(ray, its))
                return totalRadiance;

            for( int i = 0;; i++ )
            {
                /* Find the surface that is visible in the requested direction */

                if( its.mesh->isEmitter() )
                {
                    // If the last interaction was specular we add the emittance of the light source
                    if ( its.toLocal(-ray.d).z() > 0 )
                    {
                        totalRadiance += its.mesh->getEmitter()->eval({}) * beta * wBRDF;
                    }
                }

                // Direct lighting
                if( its.mesh->getBSDF()->isDiffuse() )
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
                        float distToLightSquared = vecToLight.squaredNorm();
                        float geometricTerm = fmax(0.0f, cosI * cosO) / distToLightSquared;
                        BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(dirToLight), ESolidAngle);
                        float pdfMats = its.mesh->getBSDF()->pdf(bRec);

                        // MIS calculations
                        float pdfEms = cosO > 0.0f ?
                                  emitterPdf * emitter->getEmitter()->pdfPos(eqr) * distToLightSquared / cosO
                                : 0.0f;

                        if (pdfEms < 0.0f)
                        {
                            std::cout << "pdfEms is negative" << pdfEms << std::endl;
                            return totalRadiance;
                        }


                        if (pdfMats + pdfEms == 0.0f)
                        {
                            return totalRadiance;
                        }
                        float WLight = pdfEms/(pdfEms + pdfMats);

                        totalRadiance += its.mesh->getBSDF()->eval(bRec) *
                                         beta *
                                         geometricTerm *
                                         eqr.radiance *
                                         WLight /
                                       ( emitterPdf * lightLocPdf );
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
                BSDFQueryRecord bRec(its.toLocal(-ray.d));

                beta *= its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                float pdfMats = its.mesh->getBSDF()->pdf(bRec);
                eta /= bRec.eta;
                ray = Ray3f(its.p, its.toWorld(bRec.wo));

                // Next intersection
                if(!scene->rayIntersect(ray, its))
                {
                    return totalRadiance;
                }

                float distSquared = (its.p - ray.o).squaredNorm();
//                float cosI = its.shFrame.n.dot(its.toLocal(-ray.d));
                float cosI = its.shFrame.n.dot(-ray.d);
                if( its.mesh->isEmitter() && !lastWasSpecular)
                {
                    float pdfEms = cosI > 0.0f ? scene->emitterSamplePdf() * its.mesh->getEmitter()->pdfPos({}) * distSquared / cosI : 0.0f;
                    if( pdfMats + pdfEms == 0.0f )
                    {
                        return totalRadiance;
                    }
                    wBRDF = pdfMats / (pdfMats + pdfEms);
                }
                if( lastWasSpecular )
                {
                    wBRDF = 1.0f;
                }
            }
            return totalRadiance;
        }

    public:

        MisIntigrator(const PropertyList &props) {
        }

        std::string toString() const override {
            return "MisIntigrator[]";
        }

        std::string getName() const {
            return "MIS";
        }

    };

    NORI_REGISTER_CLASS(MisIntigrator, "path_mis")
NORI_NAMESPACE_END