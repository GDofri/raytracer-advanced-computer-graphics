//
// Created by dofri on 15.3.2024.
//


#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>
NORI_NAMESPACE_BEGIN

    class AmbientOcclusionIntegrator : public SampleIntegrator {
    public:
        AmbientOcclusionIntegrator(const PropertyList &props) {
            /* No parameters this time */
        }

        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
            /* Find the surface that is visible in the requested direction */

            Intersection its;
            if (!scene->rayIntersect(ray, its))
                return Color3f(0.0f);

            // Check if light source is occluded from the intersection point
            Vector3f direction = its.shFrame.toWorld(Warp::squareToCosineHemisphere(sampler->next2D())) ;

            Ray3f shadowRay = Ray3f(its.p, direction);
            if(scene->rayIntersect(shadowRay))
            {
                return Color3f(0.0f);
            }

            // Light is not occluded.
            return Color3f(1.0);
        }

        std::string toString() const {
            return "AmbientOcclusionIntegrator[]";
        }
    };

    NORI_REGISTER_CLASS(AmbientOcclusionIntegrator, "ao");
NORI_NAMESPACE_END