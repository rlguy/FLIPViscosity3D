"""
This script will import a sequence of .obj or .ply files into Blender and
duplicate a particle object over the vertices for rendering.

Usage: copy this script into the Blender text editor and press 'Run Script'
"""

import bpy, os

# Change this to the directory containing the mesh files
MESH_DIRECTORY = "C:/example/path/to/mesh/files"

# Mesh format (either 'OBJ' or 'PLY')
MESH_FORMAT = 'OBJ'

# Mesh translation, rotation (euler, degrees), and scale 
MESH_TRANSLATION = (0, 0, 0)
MESH_ROTATION    = (0, 0, 0)
MESH_SCALE       = (10, 10, 10)

# Size of the particle object
PARTICLE_SCALE = (0.05, 0.05, 0.05)

# Internal global variables
_LOADED_MESH_NAME = ""
_LOADED_DUPLIVERT_NAME = ""
_LOADED_FRAME = -1


def remove_mesh(mesh_name):
    bl_object = bpy.data.objects.get(mesh_name)
    if bl_object is None:
        return
    mesh_data = bl_object.data
    bpy.data.objects.remove(bl_object, True)
    mesh_data.user_clear()
    bpy.data.meshes.remove(mesh_data)


def generate_icosphere(object_name):
    # Icosphere with 1 subdivision centred at origin
    verts = [
        (0.0000, 0.0000, -1.0000),
        (0.7236, -0.5257, -0.4472),
        (-0.2764, -0.8506, -0.4472),
        (-0.8944, 0.0000, -0.4472),
        (-0.2764, 0.8506, -0.4472),
        (0.7236, 0.5257, -0.4472),
        (0.2764, -0.8506, 0.4472),
        (-0.7236, -0.5257, 0.4472),
        (-0.7236, 0.5257, 0.4472),
        (0.2764, 0.8506, 0.4472),
        (0.8944, 0.0000, 0.4472),
        (0.0000, 0.0000, 1.0000)
    ]

    faces = [
        (0, 1, 2),
        (1, 0, 5),
        (0, 2, 3),
        (0, 3, 4),
        (0, 4, 5),
        (1, 5, 10),
        (2, 1, 6),
        (3, 2, 7),
        (4, 3, 8),
        (5, 4, 9),
        (1, 10, 6),
        (2, 6, 7),
        (3, 7, 8),
        (4, 8, 9),
        (5, 9, 10),
        (6, 10, 11),
        (7, 6, 11),
        (8, 7, 11),
        (9, 8, 11),
        (10, 9, 11)
    ]

    mesh = bpy.data.meshes.new(object_name + "_mesh")
    mesh.from_pydata(verts, [], faces)
    for p in mesh.polygons:
        p.use_smooth = True
    obj = bpy.data.objects.new(object_name, mesh)
    bpy.context.scene.objects.link(obj) 
    return obj


def load_mesh(scene, mesh_filepath):
    global MESH_FORMAT
    global MESH_TRANSLATION
    global MESH_ROTATION
    global MESH_SCALE
    global PARTICLE_SCALE
    global _LOADED_MESH_NAME
    global _LOADED_DUPLIVERT_NAME
    global _LOADED_FRAME

    if scene.frame_current == _LOADED_FRAME:
        return

    if _LOADED_MESH_NAME:
        remove_mesh(_LOADED_MESH_NAME)
    if _LOADED_DUPLIVERT_NAME:
        remove_mesh(_LOADED_DUPLIVERT_NAME)

    if MESH_FORMAT == 'OBJ':
        bpy.ops.import_scene.obj(filepath=mesh_filepath)
    elif MESH_FORMAT == 'PLY':
        bpy.ops.import_mesh.ply(filepath=mesh_filepath)

    mesh_object = bpy.context.selected_objects[0]
    mesh_object.location = MESH_TRANSLATION
    mesh_object.rotation_euler = MESH_ROTATION
    mesh_object.scale = MESH_SCALE

    particle_object = generate_icosphere("particle_object")
    particle_object.scale[0] *= PARTICLE_SCALE[0] * (1.0 / MESH_SCALE[0])
    particle_object.scale[1] *= PARTICLE_SCALE[1] * (1.0 / MESH_SCALE[1])
    particle_object.scale[2] *= PARTICLE_SCALE[2] * (1.0 / MESH_SCALE[2])
    particle_object.parent = mesh_object
    mesh_object.dupli_type = 'VERTS'

    _LOADED_MESH_NAME = mesh_object.name
    _LOADED_DUPLIVERT_NAME = particle_object.name
    _LOADED_FRAME = scene.frame_current


def frame_change_pre(scene):
    global MESH_DIRECTORY
    global MESH_FORMAT

    frameno = scene.frame_current
    framestr = str(frameno).zfill(4)
    if MESH_FORMAT == 'OBJ':
        filename = framestr + ".obj"
    elif MESH_FORMAT == 'PLY':
        filename = framestr + ".ply"

    filepath = os.path.join(os.path.normpath(MESH_DIRECTORY), filename)
    if not os.path.isfile(filepath):
        return
    load_mesh(scene, filepath)


bpy.app.handlers.frame_change_pre.append(frame_change_pre)