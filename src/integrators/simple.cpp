//
// Created by dofri on 15.3.2024.
//

#include <nori/integrator.h>
#include <nori/scene.h>
#include <math.h>
NORI_NAMESPACE_BEGIN

    /** SimpleIntegrator, uses simple point light source.
     */
    class SimpleIntegrator : public SampleIntegrator {
    public:
        SimpleIntegrator(const PropertyList &props) {
            lightEnergy = props.getColor("energy");
            lightPosition = props.getPoint("position");
        }

        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
            /* Find the surface that is visible in the requested direction */

            Intersection its;
            if (!scene->rayIntersect(ray, its))
                return Color3f(0.0f);

            // Check if light source is occluded from the intersection point.
            Ray3f shadowRay = Ray3f(its.p, lightPosition - its.p);
            if(scene->rayIntersect(shadowRay))
            {
                // Light is occluded
                return Color3f(0.0f);
            }

            // Light is not occluded
            Vector3f dirToLight = (lightPosition-its.p).normalized();
            float incidentCos = its.shFrame.n.dot(dirToLight);
            Vector3f vecToLight = its.p - lightPosition;
            float distToLightSquared = vecToLight.dot(vecToLight);

            Color3f color = Color3f(1.0f) = lightEnergy*fmax(0.0f, incidentCos)/(4*M_PI*M_PI * distToLightSquared);
            return color;
        }

        std::string toString() const {
            return "SimpleIntegrator[]";
        }

        std::string getName() const {
            return "Simple";
        }

        private:
        Color3f lightEnergy;
        Point3f lightPosition;
    };

    NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END