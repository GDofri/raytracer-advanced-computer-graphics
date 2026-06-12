/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/camera.h>
#include <nori/rfilter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>
NORI_NAMESPACE_BEGIN

/**
 * \brief Perspective camera with depth of field
 *
 * This class implements a simple perspective camera model. It uses an
 * infinitesimally small aperture, creating an infinite depth of field.
 */
class PerspectiveCamera : public Camera {
public:
    PerspectiveCamera(const PropertyList &propList) {
        /* Width and height in pixels. Default: 720p */
        m_outputSize.x() = propList.getInteger("width", 1280);
        m_outputSize.y() = propList.getInteger("height", 720);
        m_invOutputSize = m_outputSize.cast<float>().cwiseInverse();

        /* Specifies an optional camera-to-world transformation. Default: none */
        m_cameraToWorld = propList.getTransform("toWorld", Transform());

        /* Horizontal field of view in degrees */
        m_fov = propList.getFloat("fov", 30.0f);

        /* Near and far clipping planes in world-space units */
        m_nearClip = propList.getFloat("nearClip", 1e-4f);
        m_farClip = propList.getFloat("farClip", 1e4f);

        m_rfilter = NULL;
    }

    void activate() {
        float aspect = m_outputSize.x() / (float) m_outputSize.y();

        /* Project vectors in camera space onto a plane at z=1:
         *
         *  xProj = cot * x / z
         *  yProj = cot * y / z
         *  zProj = (far * (z - near)) / (z * (far-near))
         *  The cotangent factor ensures that the field of view is 
         *  mapped to the interval [-1, 1].
         */
        float recip = 1.0f / (m_farClip - m_nearClip),
              cot = 1.0f / std::tan(degToRad(m_fov / 2.0f));

        // The warping of z sets it to the interval [0, 1] after accounting for the homogenization (division by w)
        Eigen::Matrix4f perspective;
        perspective <<
            cot, 0,   0,   0,
            0, cot,   0,   0,
            0,   0,   m_farClip * recip, -m_nearClip * m_farClip * recip,
            0,   0,   1,   0;

        /**
         * Translation and scaling to shift the clip coordinates into the
         * range from zero to one. Also takes the aspect ratio into account.
         */
        m_sampleToCamera = Transform(
            Eigen::DiagonalMatrix<float, 3>(Vector3f(-0.5f, -0.5f * aspect, 1.0f)) *
            Eigen::Translation<float, 3>(-1.0f, -1.0f/aspect, 0.0f) * perspective).inverse();

        /* If no reconstruction filter was assigned, instantiate a Gaussian filter */
        if (!m_rfilter)
            m_rfilter = static_cast<ReconstructionFilter *>(
                NoriObjectFactory::createInstance("gaussian", PropertyList()));

        Point3f pMin = m_sampleToCamera * Point3f(0, 0, 0);
        Point3f pMax = m_sampleToCamera * Point3f(1, 1, 0);
        pMin /= pMin.z();
        pMax /= pMax.z();
        m_filmArea = std::abs((pMax.x() - pMin.x()) * (pMax.y() - pMin.y()));
        std::cout << "Camera film area: " << m_filmArea << std::endl;
        m_cameraNormal = m_cameraToWorld * Vector3f(0,0,1);
        m_cameraWorldLoc = m_cameraToWorld * Point3f(0,0,0);
        std::cout << "CameraWorldLoc: " << m_cameraWorldLoc.toString() << std::endl;
    }

    // Ray in world space.
    Color3f We( const Ray3f &ray, Point2f &rasterLoc ) const {
        float cos = std::max(0.f, ray.d.dot(m_cameraNormal));
        if(cos == 0) { return 0.f; }
        Point3f arbWorldLoc = ray(1.f);
        // Position on the image plane normalized to [0, 1]^2
        Point3f projPos = m_sampleToCamera.inverse() * m_cameraToWorld.inverse() * arbWorldLoc;
        if(projPos.x() < 0 || projPos.x() >= 1 || projPos.y() < 0 || projPos.y() >= 1)
        {
            return 0.f;
        }
        rasterLoc = Point2f(
                projPos.x() * m_outputSize.x() + 0.5f,
                projPos.y() * m_outputSize.y() + 0.5f
        );
        return 1/(m_filmArea * cos * cos * cos * cos);
    }

    void pdf(const Ray3f &ray, float &pdfPos, float &pdfDir) const override{
        float cos = std::max(0.f, ray.d.dot(m_cameraNormal));
        if(cos == 0) { pdfPos = 0.f; pdfDir = 0.f; return; }
        Point3f arbWorldLoc = ray(1.f);
        // Position on the image plane normalized to [0, 1]^2
        Point3f projPos = m_sampleToCamera.inverse() * m_cameraToWorld.inverse() * arbWorldLoc;
        if(projPos.x() < 0 || projPos.x() >= 1 || projPos.y() < 0 || projPos.y() >= 1)
        {
            pdfPos = 0.f;
            pdfDir = 0.f;
            return;
        }
        pdfPos = 1.f;
        pdfDir = 1.f / m_filmArea * (cos * cos * cos);
    }

    Color3f sampleWi(const nori::Intersection &its, Point2f &rasterLoc, float &pdf) const override {

        Point3f worldRefLoc = its.p;
        Point3f cameraRefLoc = m_cameraToWorld.inverse() * worldRefLoc;

        Point3f projPos = m_sampleToCamera.inverse() * cameraRefLoc;

        /* Check if the point lies within the image plane */
        if (projPos.x() < 0 || projPos.x() >= 1 || projPos.y() < 0 || projPos.y() >= 1)
            return 0.f;

        // Convert from [0, 1] to raster coordinates
        rasterLoc = Point2f(
                projPos.x() * m_outputSize.x(),
                projPos.y() * m_outputSize.y()
        );

        Vector3f vecToCam = (m_cameraWorldLoc - worldRefLoc);
        float distToRasterSquared = vecToCam.squaredNorm();
        float cos = std::max(0.f, m_cameraNormal.dot(-vecToCam.normalized()));
//        float cos = std::abs( m_cameraNormal.dot(-vecToCam.normalized()));
        if(cos == 0)
        {
            pdf = 0.f;
            return 0.f;
        }
        pdf = distToRasterSquared / cos;

        Ray3f ray = Ray3f(m_cameraWorldLoc, -vecToCam.normalized());
        // Importance
        return We(ray, rasterLoc);
    }

    Color3f sampleRay(Ray3f &ray,
            const Point2f &samplePosition,
            const Point2f &apertureSample) const {
        /* Compute the corresponding position on the 
           near plane (in local camera space) */
        Point3f nearP = m_sampleToCamera * Point3f(
            samplePosition.x() * m_invOutputSize.x(),
            samplePosition.y() * m_invOutputSize.y(), 0.0f);

        /* Turn into a normalized ray direction, and
           adjust the ray interval accordingly */
        Vector3f d = nearP.normalized();
        float invZ = 1.0f / d.z();

        ray.o = m_cameraToWorld * Point3f(0, 0, 0);
        ray.d = m_cameraToWorld * d;
        ray.mint = m_nearClip * invZ;
        ray.maxt = m_farClip * invZ;
        ray.update();

        return Color3f(1.0f);
    }

    void addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case EReconstructionFilter:
                if (m_rfilter)
                    throw NoriException("Camera: tried to register multiple reconstruction filters!");
                m_rfilter = static_cast<ReconstructionFilter *>(obj);
                break;

            default:
                throw NoriException("Camera::addChild(<%s>) is not supported!",
                    classTypeName(obj->getClassType()));
        }
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "PerspectiveCamera[\n"
            "  cameraToWorld = %s,\n"
            "  outputSize = %s,\n"
            "  fov = %f,\n"
            "  clip = [%f, %f],\n"
            "  rfilter = %s\n"
            "]",
            indent(m_cameraToWorld.toString(), 18),
            m_outputSize.toString(),
            m_fov,
            m_nearClip,
            m_farClip,
            indent(m_rfilter->toString())
        );
    }
private:
    Vector2f m_invOutputSize;
    Transform m_sampleToCamera;
    Transform m_cameraToWorld;
    float m_fov;
    float m_nearClip;
    float m_farClip;
    float m_filmArea;
};

NORI_REGISTER_CLASS(PerspectiveCamera, "perspective");
NORI_NAMESPACE_END
