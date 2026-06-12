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
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));
        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {

        if (  //  Frame::cosTheta(bRec.wi) <= 0 ||
        Frame::cosTheta(bRec.wo) <= 0)
            return Color3f(0.0f);

        Vector3f wh = (bRec.wi + bRec.wo).normalized();

        auto chi = [](float c){return c > 0 ? true : false;};

        auto g1 = [ this, &chi ](const Vector3f& wv, const Vector3f& wh){
            Vector3f n =  Vector3f(0,0,1);
            if( chi(wv.dot(wh)/wv.dot(n) ) )
            {
                float b = 1/(this->m_alpha*Frame::tanTheta(wv));
                return b < 1.6 ?
                       ( 3.535f*b + 2.181f*b*b ) / ( 1.0f + 2.276f*b + 2.577f*b*b) :
                       1.0f;
            }
            else
            {
                return 0.0f;
            }
        };

        float beckman = Warp::squareToBeckmannPdf(wh, m_alpha);
        float fre = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);

        float g = g1(bRec.wi, wh) * g1(bRec.wo, wh);
        float denominator = 4
                * Frame::cosTheta(bRec.wi)
                * Frame::cosTheta(bRec.wo)
                * Frame::cosTheta(wh);
        if( denominator == 0.0f)
        {
            std::cout << "dividing by zeroo" << std::endl;
        }
        return m_kd/M_PI + m_ks * beckman * fre * g / denominator;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {

        if( bRec.wo.z() <= 0.f) { return 0.0f; }
        if( bRec.wi.z() <= 0.f) { return 0.0f; }

        Vector3f wh = (bRec.wi + bRec.wo).normalized();
        float out = m_ks
            * Warp::squareToBeckmannPdf(wh, m_alpha)
            * 1 / ( 4 * wh.dot(bRec.wo) )
            + (1-m_ks) * Frame::cosTheta(bRec.wo)/M_PI;
        return out;
    }



    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {

        if( bRec.wi.z() <= 0.f) { return 0.0f; }

        if( _sample.x() <= m_ks)
        {
            // Specular case
            float reuse = _sample.x()/m_ks;
            Vector3f n = Warp::squareToBeckmann(Point2f(reuse, _sample.y()), m_alpha);
            bRec.wo = (2*n.dot(bRec.wi)*n - bRec.wi).normalized();
        }
        else
        {
            // Diffuse case
            float reuse = _sample.x()/(1-m_ks) - m_ks/(1-m_ks);
            bRec.wo = Warp::squareToCosineHemisphere(Point2f(reuse, _sample.y()));
        }
        if( bRec.wo.z() <= 0.f) { return 0.0f; }

        // Note: Once you have implemented the part that computes the scattered
        // direction, the last part of this function should simply return the
        // BRDF value divided by the solid angle density and multiplied by the
        // cosine factor from the reflection equation, i.e.
        Color3f evaled = eval(bRec);
        float cos = Frame::cosTheta(bRec.wo);
        float pdf = this->pdf(bRec);
        return evaled * cos / pdf;
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:

    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
