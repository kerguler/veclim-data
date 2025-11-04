import importlib, pkgutil

modules = {}
for _, name, _ in pkgutil.iter_modules(__path__):
    modules[name] = importlib.import_module(f"{__name__}.{name}")

tile_dat = {}
for tile in modules:
    for dat in modules[tile]:
        if dat in tile_dat:
            print("WARNING: Replacing %s with version %s!" %(dat,tile))
        tile_dat[dat] = modules[tile][dat]