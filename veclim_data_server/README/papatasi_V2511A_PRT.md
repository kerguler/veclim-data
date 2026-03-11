# Sand Fly Model API – `papatasi_V2511A_PRT`

This document describes how to access simulation outputs from the **ISMED-CLIM zoonotic disease risk model** `papatasi_V2511A_PRT` through the **VEClim data server API**.

The API returns **JSON-formatted outputs** including simulated vector population dynamics, derived risk indicators, and surveillance observations (if present) for a specified location and time period. The spatial domain is restricted to Portugal, and the temporal domain consists of climatological daily averages indexed by day-of-year values (0–364).

API endpoint:

```
https://veclim.com/api
```

---

# Example Request

A typical request can be made using `curl`:

```bash
curl "https://veclim.com/api?vec=papatasi_V2511A_PRT&lat=38.25&lon=-8.25&dates=2026-03-11:2026-04-11&opr=ts"
```

---

# Request Parameters

### `vec`

Name of the vector model.

For the sand fly model described here:

```
vec=papatasi_V2511A_PRT
```

---

### `lat`

Latitude of the requested location (decimal degrees).

Example:

```
lat=38.25
```

---

### `lon`

Longitude of the requested location (decimal degrees).

Example:

```
lon=-8.25
```

---

### `date` or `dates`

Defines the requested time period.

Accepted formats:

Single day

```
date=YYYY-MM-DD
```

Example:

```
date=2026-03-11
```

Time interval

```
dates=YYYY-MM-DD:YYYY-MM-DD
```

Example:

```
dates=2026-03-11:2026-04-11
```

---

### `opr`

Defines the type of returned data.

Accepted values:

| value  | description                                   |
| ------ | --------------------------------------------- |
| `ts`   | daily time series                             |
| `mean` | daily climatological mean or temporal average (not applicable for this model) |

Example:

```
opr=ts
```

---

# Response Structure

The API response is a JSON object containing several sections.

## `location`

Information about the spatial unit corresponding to the requested coordinates. The spatial units are in reference to [**georef-portugal-concelho-millesime**](https://public.opendatasoft.com/explore/assets/georef-portugal-concelho-millesime/).

Fields may include:

```
lon      longitude
lat      latitude
pid      polygon identifier
name     location name
island   island indicator
```

Example:

```json
"location": {
  "lon": -8.25,
  "lat": 38.25,
  "pid": "1501",
  "name": "alcácer do sal",
  "island": 1
}
```

---

## `date`

Information about the requested time period.

```
date0   start date
date1   end date
days    day-of-year interval (0..364)
valid   request validity indicator
```

---

## `sim-ts`

Time series of **model simulation outputs**.

For model `V2511A_PRT` the following variable is currently provided:

| variable | description                                   |
| -------- | --------------------------------------------- |
| `newegg` | simulated daily number of newly produced eggs as a proxy to sand fly activity |

Example:

```json
"sim-ts": {
  "V2511A_PRT": {
    "newegg": [ ... ]
  }
}
```

---

## `risk-ts`

Derived **seasonal risk indicators** computed from the simulated population dynamics.

| field       | description                              |
| ----------- | ---------------------------------------- |
| `up`        | start of the active season (day-of-year) |
| `down`      | end of the active season (day-of-year)   |
| `peak_up`   | start of peak activity period            |
| `peak_down` | end of peak activity period              |
| `risk`      | daily risk indicator (0, 1, 2)           |

Example:

```json
"risk-ts": {
  "V2511A_PRT": {
    "up": [117],
    "down": [332],
    "peak_up": [189],
    "peak_down": [290],
    "risk": [ ... ]
  }
}
```

---

## `surv-ts`

Optional surveillance data associated with the requested location.

| field        | description                             |
| ------------ | --------------------------------------- |
| `adult_norm` | normalized adult abundance observations (per sampling per day) |

Example:

```json
"surv-ts": {
  "adult_norm": []
}
```

---

# Notes

* Coordinates are automatically mapped to the nearest spatial unit in the VEClim database.
* Time series are returned at **daily resolution**.
* Seasonal indicators (`up`, `down`, `peak_up`, `peak_down`) are expressed as **day-of-year** values (0..364).
