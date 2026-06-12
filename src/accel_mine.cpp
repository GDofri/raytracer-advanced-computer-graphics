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

#include <nori/accel.h>
#include <Eigen/Geometry>
#include <chrono>
NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {

    if(m_meshes.empty())
    {
        m_bbox = mesh->getBoundingBox();
    }
    else
    {
        m_bbox.expandBy(mesh->getBoundingBox());
    }
    m_meshes.push_back(mesh);
    m_bboxes.push_back(mesh->getBoundingBox());
}

void Accel::build() {
    /* Nothing to do here for now */

    for( size_t i = 0; i < m_meshes.size(); ++i)
    {
        std::cout << "Start building octtree no. " << i << std::endl;
        std::cout << m_meshes[i]->getName() << std::endl;
        m_octtrees.push_back(std::make_unique<Octtree>(m_bboxes[i], *m_meshes[i]));
        std::cout << "Built octtree no. " << i << std::endl;
    }

//    const auto &[interior, leaf, triangles] = m_octtree->getTreeStatistics();
////    std::cout << "Tree statistics:" << std::endl;
////    std::cout << "Number of interior nodes: " << interior << std::endl;
////    std::cout << "Number of leaf nodes: " << leaf << std::endl;
////    std::cout << "Avg. number of triangles per leaf node: " << (double)triangles/leaf << std::endl;

}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {

    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    // Find the closest intersection point
    float closestT = std::numeric_limits<float>::max();
    std::optional<Intersection> closestIntersection = std::nullopt;
    uint32_t closestIndex = -1;
    // Loop through all octtrees and find the closest intersection
    // TODO Maybe later it would be better to keep all the meshes in one octtree
    for( const auto &octtree : m_octtrees  )
    {
        const auto &[index, intersection] = octtree->findClosestIntersection(ray, shadowRay);
        if( !intersection.has_value() ) continue;
        if( intersection->t < closestT )
        {
            closestT = intersection->t;
            closestIntersection = intersection;
            closestIndex = index;
            ray.maxt = closestT;

            // If it is a shadow ray we don't need to search for more than one intersection.
            if( shadowRay ) break;
        }
    }
    foundIntersection = closestIntersection.has_value();

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */



        its = closestIntersection.value();
        f = closestIndex;

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                     bary.y() * UV.col(idx1) +
                     bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                    (bary.x() * N.col(idx0) +
                     bary.y() * N.col(idx1) +
                     bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

NORI_NAMESPACE_END

