# daedong-detailmap

[한국어](#한국어) | [English](#english)

---

## 한국어

수묵화 느낌의 한반도 지도 앱 「대동명소도」에서 지도를 크게 키우면 나타나는 상세 지도 자료를 만드는 도구 모음입니다.
OpenStreetMap과 Copernicus GLO-30 높낮이 자료를 읽어서 앱이 쓰는 이진 파일을 만듭니다.
앱 코드, 원본 자료, 만들어진 자료는 이 저장소에 들어 있지 않습니다. 자료를 만드는 방법만 공개합니다.

### 만들어 내는 것

모두 한국 땅 쪽으로만 잘라 냅니다. 북쪽 경계는 압록강과 두만강의 중심선입니다.

| 내용 | 파일 |
| --- | --- |
| 배율 띠별 봉우리 산 그림 목록 (A 약 4배, B 약 8배, C 약 16배) | `out/bundle/peaks_A.bin`, `peaks_B.bin`, `peaks_C.bin` |
| 해안선과 강 (A), 개울 (B), 호수와 넓은 강 | `out/bundle/lines_A.bin`, `lines_B.bin` |
| 산줄기를 잇는 희미한 번짐용 큰 분수령 | `out/bundle/ridges.bin` |
| 한국 땅 면 (조각별) | `out/bundle/land_A.bin`, `land_B.bin` |
| 상세 층이 그릴 구역 (한국 땅과 바다 25 km 안) | `out/bundle/area.bin` |
| 이름표: 봉우리, 강, 개울, 시·군·읍 이름 (한국어) | `out/bundle/labels.bin` |
| 이름표에 필요한 글자, 글꼴에 없는 글자 | `out/bundle/chars.txt`, `chars_missing.txt` |
| 일반 지도의 손그림 산을 실제 봉우리로 옮기는 값 | `out/general/glyph_moves.json`, `MapBase_moved.json` |
| 일반 지도용 한국 땅 해안 (OSM 해안선) | `out/general/korea_land_general.json` |
| 검사 값 (해안선 차이, 파일 크기) | `out/land_check.json`, `out/sizes.json` |

파일 짜임새는 `bundle.py` 맨 위 설명과 `geom.py`의 `pack`에 적혀 있습니다.
좌표는 모두 지도 점(1360 x 2400, 메르카토르 투영, 경도 124.0~131.95, 위도 32.85~43.15)입니다. 위경도에서 바꾸는 식은 `common.py`의 `proj`입니다.

### 준비물

- Python 3.9 이상, zsh (`build.sh`가 zsh 스크립트입니다. macOS는 기본으로 있습니다.)
- 파이썬 꾸러미: `requirements.txt` (numpy, rasterio, scipy, Pillow)
- 디스크 여유 약 2.5 GB (내려받는 자료 약 1.7 GB, 중간 결과 `work/` 약 200 MB, 결과 `out/` 약 10 MB)

```sh
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt
```

### 내려받을 자료

이 저장소에는 없습니다. 모두 `data/` 아래에 둡니다.

1. **Copernicus GLO-30 높낮이 조각 50장**
   `./build.sh tiles`가 내려받을 주소(AWS의 `copernicus-dem-30m`, 로그인 없이 받을 수 있음)를 한 줄씩 보여 줍니다.
   파일 이름 그대로 `data/dem/`에 둡니다.
   ```sh
   mkdir -p data/dem data/osm
   (cd data/dem && ../../build.sh tiles | xargs -n 1 -P 4 curl -fLO)
   ```
   자료 소개: https://registry.opendata.aws/copernicus-dem/
2. **Geofabrik의 대한민국, 북한 PBF** (`data/osm/`)
   ```sh
   curl -fL -o data/osm/south-korea-latest.osm.pbf https://download.geofabrik.de/asia/south-korea-latest.osm.pbf
   curl -fL -o data/osm/north-korea-latest.osm.pbf https://download.geofabrik.de/asia/north-korea-latest.osm.pbf
   ```
3. **일반 지도 장면** (아래 "일반 지도 장면" 참고)

### 실행

```sh
export DETAILMAP_MAPBASE=/경로/MapBase.json   # 또는 data/MapBase.json 에 둔다
./build.sh
```

`build.sh`는 일반 지도 장면이 없으면 처음에 멈추고 이유를 알려 줍니다. 단계는 다음과 같습니다.

| 단계 | 하는 일 |
| --- | --- |
| `osm_extract.py` | PBF에서 해안, 물줄기, 호수, 이름을 뽑습니다. (`osmpbf.py`는 추가 꾸러미 없는 PBF 읽개) |
| `osm_boundary.py` | 남북의 영역을 뽑습니다. 국경이 강을 따르는 곳의 땅 경계가 됩니다. |
| `korea_mask.py` | 일반 지도 장면으로 한국 땅 가리개를 만듭니다. |
| `terrain.py` | 높낮이 조각마다 두드러진 봉우리와 큰 분수령을 찾습니다. (조각마다 따로 돌며 5개씩 동시에 실행) |
| `land.py` | 한국 땅 면, 자르기 가리개, 그릴 구역, 일반 지도용 해안을 만듭니다. |
| `bundle.py` | 위의 결과를 묶어 `out/bundle/`의 이진 파일과 이름표, 글자 목록을 씁니다. |

중간 결과는 `work/`, 최종 결과는 `out/`에 쌓입니다. `./build.sh tiles`만 쓰면 내려받을 주소만 보여 주고 끝납니다.

원래 저장소에는 앱의 그리기 코드로 결과를 눈으로 확인하는 단계(`preview`)가 더 있었지만, 앱 코드라서 이 저장소에는 넣지 않았습니다.
`build.sh`는 그 파일이 없으면 이 단계를 알아서 건너뜁니다.

### 일반 지도 장면 (이 저장소에 없음)

앱의 일반 지도(한반도 전체)는 JSON 장면 파일 `MapBase.json`으로 그려집니다. 이 파일은 포함하지 않았습니다.
`korea_mask.py`, `land.py`, `bundle.py`가 이 파일을 읽으므로, 위치를 환경변수 `DETAILMAP_MAPBASE`로 알려 주거나 `data/MapBase.json`에 두어야 합니다.
없으면 `build.sh`와 각 스크립트가 분명한 오류로 멈춥니다. (`korea_mask.py`는 `python -I korea_mask.py <MapBase.json> <out.npz>` 꼴로 경로를 직접 줄 수도 있습니다.)

같은 짜임새라면 어떤 장면으로 바꿔 써도 됩니다. 장면에서 실제로 쓰는 것은 아래뿐입니다.

| 장면 안의 것 | 쓰는 곳 | 쓰임 |
| --- | --- | --- |
| `geoms.korea` : 한국 땅의 닫힌 고리 목록. 고리 하나는 `[x0, y0, x1, y1, ...]` 숫자 목록이고, 좌표는 위의 지도 점입니다. | `korea_mask.py`, `land.py` | 가리개를 만듭니다. (고리 안쪽을 채우고 지도 점 1.5개, 약 0.8 km 넓힘.) 가리개는 압록강, 두만강 건너 중국 쪽 봉우리와 산줄기를 버리는 데 씁니다. `land.py`는 이 해안이 OSM 해안과 얼마나 다른지 `out/land_check.json`에 적습니다. |
| `geoms.korea2` : 한국과 땅으로 이어진 이웃 땅의 고리 목록 (같은 형식) | `land.py` | 가장 긴 고리의 국경 쪽 호를 써서 일반 지도용 해안을 닫습니다. (`korea_land_general.json`) |
| `ops` : 그리기 명령 목록. 그중 처음으로 나오는 `t`가 `"glyphs"`이고 `items`가 비어 있지 않은 명령. 항목은 `p`(숫자 14개 = 점 7개)와 선택 항목 `x`, `d`를 가집니다. | `bundle.py` | 일반 지도의 손그림 산을 가까운 실제 봉우리로 옮깁니다. |
| `detailApplied` : 참이면 이미 옮긴 장면으로 보고, `out/general/glyph_moves.json`이 이미 있어야 합니다. | `bundle.py` | 없어도 됩니다. |

그래서 일반 지도 장면이 따로 없다면 다음처럼 직접 만든 가짜 장면으로도 돌아갑니다.

- 해안선 고리를 `common.proj`로 지도 점으로 바꿔 `geoms.korea`와 `geoms.korea2`에 넣고,
- `ops`에는 산 그림 항목이 하나뿐인 명령을 넣습니다.
  `{"t": "glyphs", "items": [{"p": [100,100, 101,99, 102,98, 103,97, 104,98, 105,99, 106,100]}]}`

이 경우 `out/general/`의 일반 지도용 결과(옮김 값, 해안)는 쓸모가 없지만, `out/bundle/`의 결과는 바뀌지 않습니다.

`bundle.py`는 글꼴에 없는 글자(`chars_missing.txt`)를 찾을 때 한글 완성형 2350자와 ASCII를 기준으로 삼습니다. 앱 소스(`../../Daedong`)가 옆에 있으면 거기 쓰인 글자도 기준에 더하고, 없으면 그냥 건너뜁니다.

### 사용권과 출처

- **코드**: MIT 사용권입니다. [LICENSE](LICENSE)를 보세요.
- **지도 자료**: 이 도구로 만든 자료는 OpenStreetMap에서 파생된 자료입니다.
  © OpenStreetMap contributors. [Open Database License (ODbL) 1.0](https://opendatacommons.org/licenses/odbl/1-0/)으로 쓸 수 있습니다.
  자세한 것은 https://www.openstreetmap.org/copyright 에 있습니다.
  이 저장소는 파생 자료를 만든 방법을 누구나 볼 수 있게 하는 것이기도 합니다. (ODbL 제4.6조)
- **높낮이**: 아래 문장을 그대로 밝혀야 합니다.

  > produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved.

- 일반 지도 장면은 이 저장소에 없으므로 쓰는 장면의 사용권은 따로 확인하세요.

---

## English

Tools that build the "detail map" data for the iOS app 대동명소도 (Daedong Myeongsodo), an ink-wash style map of Korea. When the user zooms in, the app draws real mountain summits as brush-style glyphs, with coastline, rivers, streams, lakes and Korean place names on top. These scripts read OpenStreetMap and Copernicus GLO-30 elevation data and write the binary files the app uses.
The app code, the source data and the generated data are not in this repository. Only the method of producing the data is published.

### What it produces

Everything is clipped to the Korean side; the northern limit is the centre line of the Yalu (압록강) and Tumen (두만강) rivers.

| Content | File |
| --- | --- |
| Peak glyph lists per zoom band (A about 4x, B about 8x, C about 16x) | `out/bundle/peaks_A.bin`, `peaks_B.bin`, `peaks_C.bin` |
| Coastline and rivers (A), streams (B), lakes and wide rivers | `out/bundle/lines_A.bin`, `lines_B.bin` |
| Major drainage divides for the faint ridge wash between glyphs | `out/bundle/ridges.bin` |
| Korean land polygons, per tile | `out/bundle/land_A.bin`, `land_B.bin` |
| Area the detail layer draws (Korean land plus sea within 25 km) | `out/bundle/area.bin` |
| Labels: peaks, rivers, streams, city/county/town names (Korean) | `out/bundle/labels.bin` |
| Characters the labels need, and those missing from the font | `out/bundle/chars.txt`, `chars_missing.txt` |
| Moves that put the general map's hand-drawn mountains on real summits | `out/general/glyph_moves.json`, `MapBase_moved.json` |
| Korean land for the general map (OSM coastline) | `out/general/korea_land_general.json` |
| Check values (coastline offset, file sizes) | `out/land_check.json`, `out/sizes.json` |

The file formats are described at the top of `bundle.py` and in `geom.pack`.
All coordinates are map points: a 1360 x 2400 plane in Mercator projection covering longitude 124.0 to 131.95 and latitude 32.85 to 43.15. `proj` in `common.py` converts longitude and latitude.

### Requirements

- Python 3.9 or newer, and zsh (`build.sh` is a zsh script; macOS has it by default)
- Python packages in `requirements.txt` (numpy, rasterio, scipy, Pillow)
- About 2.5 GB of free disk (downloads about 1.7 GB, intermediate `work/` about 200 MB, results `out/` about 10 MB)

```sh
python3 -m venv .venv
./.venv/bin/pip install -r requirements.txt
```

### Downloads

None of these are in the repository. Put them under `data/`.

1. **Copernicus GLO-30 DEM tiles, 50 files.**
   `./build.sh tiles` prints one URL per line (AWS bucket `copernicus-dem-30m`, no login needed). Keep the file names and save them in `data/dem/`.
   ```sh
   mkdir -p data/dem data/osm
   (cd data/dem && ../../build.sh tiles | xargs -n 1 -P 4 curl -fLO)
   ```
   Dataset page: https://registry.opendata.aws/copernicus-dem/
2. **Geofabrik PBF extracts for South Korea and North Korea** in `data/osm/`:
   ```sh
   curl -fL -o data/osm/south-korea-latest.osm.pbf https://download.geofabrik.de/asia/south-korea-latest.osm.pbf
   curl -fL -o data/osm/north-korea-latest.osm.pbf https://download.geofabrik.de/asia/north-korea-latest.osm.pbf
   ```
3. **A general-map scene** (see "General-map scene" below).

### Running

```sh
export DETAILMAP_MAPBASE=/path/to/MapBase.json   # or put the file at data/MapBase.json
./build.sh
```

`build.sh` stops at the start, with an explanation, if the general-map scene is missing. The steps:

| Step | What it does |
| --- | --- |
| `osm_extract.py` | Pulls coastline, waterways, lakes and names out of the PBF files (`osmpbf.py` is a PBF reader that needs no extra packages). |
| `osm_boundary.py` | Extracts the national territories of the two countries; this gives the land border where it follows the border rivers. |
| `korea_mask.py` | Builds a Korean-land mask from the general-map scene. |
| `terrain.py` | Finds prominent summits and major divides per DEM tile (one process per tile, five at a time). |
| `land.py` | Builds the Korean land polygons, the clipping masks, the draw area and the coast for the general map. |
| `bundle.py` | Packs everything into the binary files, labels and character lists in `out/bundle/`. |

Intermediate results go to `work/`, final results to `out/`. `./build.sh tiles` only prints the download URLs and exits.

The original repository has one more step (`preview`) that renders the results with the app's own drawing code. That is app code, so it is not published here. `build.sh` skips this step when its files are absent.

### General-map scene (not included)

The app's general map (the whole peninsula) is drawn from a JSON scene file, `MapBase.json`. It is not included.
`korea_mask.py`, `land.py` and `bundle.py` read it, so give its location in the environment variable `DETAILMAP_MAPBASE` or put it at `data/MapBase.json`.
If it is missing, `build.sh` and each script stop with a clear error. (`korea_mask.py` also accepts the path directly: `python -I korea_mask.py <MapBase.json> <out.npz>`.)

Any scene with the same structure can stand in. This is all that is used from it:

| In the scene | Used by | Purpose |
| --- | --- | --- |
| `geoms.korea`: list of closed rings of Korean land. Each ring is a flat `[x0, y0, x1, y1, ...]` list in the map points above. | `korea_mask.py`, `land.py` | Filled and widened by 1.5 map points (about 0.8 km) to make the mask that drops summits and ridges on the Chinese side of the Yalu and Tumen. `land.py` also writes how far this coast is from the OSM coast to `out/land_check.json`. |
| `geoms.korea2`: list of rings (same format) of the neighbouring land that joins Korea over land. | `land.py` | The border-side arc of the longest ring closes the general map's coast. (`korea_land_general.json`) |
| `ops`: the list of draw commands; the first one with `t == "glyphs"` and non-empty `items`. Each item has `p` (14 numbers = 7 points) and optional `x`, `d`. | `bundle.py` | Moves the general map's hand-drawn mountains onto nearby real summits. |
| `detailApplied`: if true, the scene is treated as already moved and `out/general/glyph_moves.json` must already exist. | `bundle.py` | Optional. |

So without the original scene a hand-made stand-in works:

- put your coastline rings, converted with `common.proj`, into `geoms.korea` and `geoms.korea2`;
- put one command with a single glyph item into `ops`:
  `{"t": "glyphs", "items": [{"p": [100,100, 101,99, 102,98, 103,97, 104,98, 105,99, 106,100]}]}`

With such a stand-in the general-map results in `out/general/` (moves, coast) are not useful, but the files in `out/bundle/` do not change.

To find characters missing from the font (`chars_missing.txt`), `bundle.py` compares against the 2,350 KS X 1001 Hangul syllables plus ASCII. If the app sources (`../../Daedong`) sit next to this folder their characters are added to that set; otherwise that part is skipped.

### Licences and attribution

- **Code**: MIT, see [LICENSE](LICENSE).
- **Map data**: data produced with these tools is derived from OpenStreetMap.
  © OpenStreetMap contributors, available under the [Open Database License (ODbL) 1.0](https://opendatacommons.org/licenses/odbl/1-0/).
  See https://www.openstreetmap.org/copyright.
  This repository also serves as the published method of making the derived data (ODbL section 4.6).
- **Elevation**: the following attribution must be given as is:

  > produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved.

- The general-map scene is not part of this repository; check the licence of whatever scene you use.
