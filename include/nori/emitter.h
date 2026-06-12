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

#pragma once

#include <nori/object.h>
#include <nori/mesh.h>
NORI_NAMESPACE_BEGIN



/**
 * \brief Convenience data structure used to pass multiple
 * parameters to the evaluation and sampling routines in \ref Emitter
 */
struct EmitterQueryRecord {

    /// Location on the surface of the emitter
    Point3f loc;

    /// Normal at the location on the surface of the emitter
    Vector3f normal;

    /// Radiance emitted by the emitter accounting for the likelihood of the sample.
    Color3f radiance;

    /// Outgoing direction from the emitter in world coordinates
    Vector3f woWorld;

    /// Probability density function of the location on the emitter
    float pdfLoc;

    /// Probability density function of the direction from the emitter
    float pdfDir;

};


/**
 * \brief Superclass of all emitters
 */
class  Emitter : public NoriObject {
public:

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }

    /**
     * \brief Set the mesh associated with this emitter
     * @param mesh
     */
    void setMesh(Mesh* mesh) { m_mesh = mesh; };

    /**
     * \brief Return the radiance emitted by the emitter
     * @return Radiance emitted by the emitter
     */
    virtual Color3f getRadiance() const = 0;

    /**
     * \brief Return the probability density function of the emitter
     *  at the point specified in the query record.
     *  @param record
     *  @return pdf
     */
    virtual float pdfPos(const EmitterQueryRecord &record) const = 0;

    /**
     * \brief Return the probability density function of the direction
     *  from the emitter in the direction specified in the query record.
     *  @param record
     *  @return pdf
     */
    virtual float pdfDir(const EmitterQueryRecord &record) const = 0;

    /**
     * \brief Evaluate the radiance emitted by the emitter
     *  at the point specified in the query record.
     *  @param record
     *  @return radiance
     */
    virtual Color3f eval(const EmitterQueryRecord &record) const = 0;

    /**
     * \brief Sample a point on the emitter
     * @param sample
     * @return EmitterQueryRecord
     */
    virtual void sample(EmitterQueryRecord &eqr, Sampler &sample) = 0;

    /**
     * Sample an outgoing ray from the emitter.
     * @param eqr
     * @param sample
     * @param ray
     */
    virtual void sampleRay(EmitterQueryRecord &eqr, Sampler &sample, Ray3f &ray) = 0;

protected:
    Mesh *m_mesh;
};

NORI_NAMESPACE_END
