import importlib, pkgutil, sys

modules = {}
for _, name, _ in pkgutil.iter_modules(__path__):
    modules[name] = importlib.import_module(f"{__name__}.{name}")

def reload_forecast_var():
    global modules
    name = "forecastECMWF"
    full = f"{__name__}.{name}"

    if full in sys.modules:
        importlib.reload(sys.modules[full])
    else:
        modules[name] = importlib.import_module(full)