"""
This script will import a sequence of .obj files and update the mesh data on a frame change.
The sequence can be rendered by converting the object to a Point Cloud using the Geometry Nodes
Mesh to Points Node.

Usage: copy this script into the Blender text editor and press 'Run Script'
"""

import bpy, os


# Change this to the directory containing the mesh files
MESH_CACHE_PATH = "C:/example/path/to/mesh/files"

_MESH_CACHE_NAME = "FLIPViscosity3D"


def frame_change_post(scene, depsgraph=None):
    original_active_object = bpy.context.view_layer.objects.active
    
    if bpy.data.objects.get(_MESH_CACHE_NAME) is None:
        mesh_data = bpy.data.meshes.new(_MESH_CACHE_NAME + "Mesh")
        mesh_data.from_pydata([], [], [])
        mesh_cache_object = bpy.data.objects.new(_MESH_CACHE_NAME, mesh_data)
        bpy.context.collection.objects.link(mesh_cache_object)
        
    mesh_cache_object = bpy.data.objects.get(_MESH_CACHE_NAME)
        
    frameno = scene.frame_current
    mesh_filename = str(frameno).zfill(4) + ".obj"
    mesh_filepath = os.path.normpath(os.path.join(MESH_CACHE_PATH, mesh_filename))
    if not os.path.isfile(mesh_filepath):
        print("File not found: <" + mesh_filepath + ">")
        return
    
    bpy.ops.wm.obj_import(filepath=mesh_filepath)
    import_object = bpy.context.view_layer.objects.active
    
    vertices = []
    for vertex in import_object.data.vertices:
        v = import_object.matrix_world @ vertex.co
        vertices.append((v[0], v[1], v[2]))
        
    import_object_mesh_data = import_object.data
    bpy.data.objects.remove(import_object, do_unlink=True)
    import_object_mesh_data.user_clear()
    bpy.data.meshes.remove(import_object_mesh_data)
        
    mesh_cache_object.data.clear_geometry()
    mesh_cache_object.data.from_pydata(vertices, [], [])
    
    bpy.context.view_layer.objects.active = original_active_object


bpy.app.handlers.frame_change_post.append(frame_change_post)