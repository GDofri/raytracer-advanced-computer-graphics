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

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {

        float probReflected = fresnel(bRec.wi.z(), m_extIOR, m_intIOR);
        bRec.measure = EDiscrete;

        if(sample.x() < probReflected)
        {
            // Reflect
            bRec.wo = Vector3f(
                -bRec.wi.x(),
                -bRec.wi.y(),
                 bRec.wi.z()
            );
            bRec.eta = 1.0f;
            return Color3f(1.0);
        }

        // Refract
        float eta = m_extIOR/m_intIOR;

        // If ray is coming to the back of a face we are "inside" the interior medium.
        Vector3f n(0.0f,0.0f,1.0f);
        if( bRec.wi.z() < 0.0f )
        {
            n *= -1.0f;
            eta = 1/eta;
        }
        Vector3f wi = bRec.wi;
        bRec.wo = -eta * ( wi -  wi.dot(n)*n) - n * sqrt( 1.0f - eta*eta * ( 1.0 - wi.dot(n)*wi.dot(n)) );
        bRec.wo.normalize();
        bRec.eta = eta;
        return Color3f(1.0);
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
