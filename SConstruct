#!/usr/bin/env python
import os
import sys

env = SConscript("godot-cpp/SConstruct")

env.Append(CPPPATH=[
    "src",
    "h3/src/h3lib/include",
    "h3/build/src/h3lib/include",
    "godot-cpp/include",
    "godot-cpp/gen/include"
])

env.Append(LIBPATH=[
    "h3/build/bin/Release"
])

env.Append(LIBS=["h3"])

env.Append(CPPPATH=["src/"])
sources = Glob("src/*.cpp")

if env["target"] in ["editor", "template_debug"]:
    try:
        doc_data = env.GodotCPPDocData("src/gen/doc_data.gen.cpp", source=Glob("doc_classes/*.xml"))
        sources.append(doc_data)
    except AttributeError:
        print("Not including class reference as we're targeting a pre-4.3 baseline.")

if env["platform"] == "macos":
    library = env.SharedLibrary(
        "h3-godot/bin/libh3_extension.{}.{}.framework/libh3_extension.{}.{}".format(
            env["platform"], env["target"], env["platform"], env["target"]
        ),
        source=sources,
    )
elif env["platform"] == "ios":
    if env["ios_simulator"]:
        library = env.StaticLibrary(
            "h3-godot/bin/libh3_extension.{}.{}.simulator.a".format(env["platform"], env["target"]),
            source=sources,
        )
    else:
        library = env.StaticLibrary(
            "h3-godot/bin/libh3_extension.{}.{}.a".format(env["platform"], env["target"]),
            source=sources,
        )
elif env["platform"] == "linux":
    env.Append(CCFLAGS=["-fPIC"])
    env.Append(LIBS=["dl"])
    library = env.SharedLibrary(
        "h3-godot/bin/libh3_extension.so",
        source=sources,
    )
elif env["platform"] == "android":
    env.Append(CCFLAGS=["-fPIC"])
    library = env.SharedLibrary(
        "h3-godot/bin/libh3_extension.so",
        source=sources,
    )
elif env["platform"] == "web":
    env.Append(CCFLAGS=["-s", "SIDE_MODULE=1"])
    env.Append(LINKFLAGS=["-s", "SIDE_MODULE=1"])
    library = env.SharedLibrary(
        "h3-godot/bin/libh3_extension.wasm",
        source=sources,
    )
else:
    library = env.SharedLibrary(
        "h3-godot/bin/libh3_extension{}{}".format(env["suffix"], env["SHLIBSUFFIX"]),
        source=sources,
    )

Default(library)
 