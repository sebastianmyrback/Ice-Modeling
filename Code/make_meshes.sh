#!/usr/bin/env bash
type gmsh >/dev/null 2>&1 || { echo >&2 "Gmsh is required to make meshes. See: http://gmsh.info/. Aborting."; exit 1; }
type meshio-convert >/dev/null 2>&1 || { echo >&2 "meshio-convert is required to make meshes. try: pip install --user meshio."; exit 1; }

version=$(gmsh --version 2>&1)
maj_version=${version:0:1}
min_version=${version:2:1}
dev_version=${version:4:1}
if [ ${maj_version} -lt 4 ] ||  [ ${min_version} -lt 4 ]; then
    echo "Gmsh version must be > 4.4.x. Your current version is ${version}"
    exit
fi

cd meshes
gmsh -3 -format msh2 extruded.geo
meshio-convert --prune extruded.msh extruded.xdmf
cd ..
