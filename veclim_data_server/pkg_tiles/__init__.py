import importlib, pkgutil

modules = {}
for _, name, _ in pkgutil.iter_modules(__path__):
    modules[name] = importlib.import_module(f"{__name__}.{name}")

tile_dat = {}
for tile in modules:
    if tile in tile_dat:
        print("WARNING: Replacing %s!" %(tile))
    tile_dat[tile] = modules[tile].tile_dat
