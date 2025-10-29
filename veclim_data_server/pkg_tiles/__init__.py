import importlib, pkgutil

modules = {}
for _, name, _ in pkgutil.iter_modules(__path__):
    modules[name] = importlib.import_module(f"{__name__}.{name}")