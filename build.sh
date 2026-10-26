#!/bin/zsh
# 대동명소도 상세 지도 자료 만들기 (detail map build)
#
# 상세 확대에서 실제 봉우리 자리에 일반 지도와 같은 붓말의 산 그림을 놓는 데 쓰는 자료를 만든다.
# 결과는 out/bundle (앱에 넣을 이진 자료), out/general (일반 지도 산 그림 옮김).
# 이 폴더의 work/ 와 out/ 밖에는 아무것도 쓰지 않는다.
#
# 준비 (한 번):
#   python3 -m venv .venv && ./.venv/bin/pip install -r requirements.txt
#   data/dem/  Copernicus GLO-30 1도 조각 50장 (AWS copernicus-dem-30m, build.sh tiles 로 목록을 본다)
#   data/osm/  Geofabrik south-korea-latest.osm.pbf, north-korea-latest.osm.pbf
# 차례:
#   ./build.sh            전부 (OSM 뽑기 -> 지형 -> 땅 -> 자료 묶기)
#   ./build.sh tiles      내려받을 Copernicus 조각 주소 목록
# 일반 지도 장면 (MapBase.json, 이 저장소에는 없음): 환경변수 DETAILMAP_MAPBASE에 경로를 주거나 data/MapBase.json에 둔다.
#   어느 것도 없으면 처음에 멈춘다. 장면에 무엇이 들어 있어야 하는지는 README.md.
# 단계:
#   osm_extract.py  PBF에서 해안, 물줄기, 호수, 이름 (osmpbf.py: 추가 패키지 없는 PBF 읽개)
#   osm_boundary.py 남북의 영역 (국경 강을 따르는 땅 경계)
#   terrain.py      조각마다 두드러짐 봉우리와 큰 분수령 (korea_mask.py로 한국 땅만)
#   land.py         한국 땅 면 (해안선 + 압록강·두만강 중심선), 자르기 가리개, 상세 층이 그릴 구역, 일반 지도용 해안
#   bundle.py       배율별 산 그림 목록, 선, 땅, 분수령 번짐, 구역, 이름표, 글자 목록, 일반 지도 산 그림 옮김, 크기
#   (확인 그림 단계 preview.py는 이 공개 판에 없다. 파일이 있을 때만 도는 칸은 남겨 두었다.)
# 결과 쓰기: out/bundle/*.bin을 앱의 자료 폴더에 복사한다. 일반 지도 쪽은 out/general/ 의 JSON을 쓴다.
# 사용권: 지도 자료는 © OpenStreetMap 기여자 (ODbL 1.0). 이 폴더의 도구가 파생 자료를 만드는 방법이다.
#   높낮이: produced using Copernicus WorldDEM-30 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH
#   2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved.
set -e
cd ${0:A:h}
PY=./.venv/bin/python
TILES=(33,126 34,125 34,126 34,127 34,128 35,125 35,126 35,127 35,128 35,129 36,125 36,126 36,127 36,128 36,129
       37,124 37,125 37,126 37,127 37,128 37,129 37,130 37,131 38,124 38,125 38,126 38,127 38,128 39,124 39,125
       39,126 39,127 39,128 40,124 40,125 40,126 40,127 40,128 40,129 41,125 41,126 41,127 41,128 41,129 41,130
       42,127 42,128 42,129 42,130 43,129)
if [[ "$1" == tiles ]]; then
  for t in $TILES; do la=${t%,*}; lo=${t#*,}; n=Copernicus_DSM_COG_10_N${la}_00_E${lo}_00_DEM
    echo https://copernicus-dem-30m.s3.amazonaws.com/$n/$n.tif; done; exit 0
fi
mkdir -p work/terrain out
HAVE_PREVIEW=0
if [[ -f preview.py && -f preview/render/main.swift && -f preview/render/MapSceneGlyph.swift ]]; then HAVE_PREVIEW=1; fi
if [[ "$1" == preview && $HAVE_PREVIEW == 0 ]]; then
  echo "preview.py 또는 preview/render가 없어 확인 그림을 만들 수 없다." >&2; exit 1
fi
if [[ $HAVE_PREVIEW == 1 ]]; then
  mkdir -p bin
  swiftc -O -o bin/render preview/render/main.swift preview/render/MapSceneGlyph.swift
fi
if [[ "$1" != preview ]]; then
  $PY -I -c 'import sys; sys.path.insert(0, "."); from common import mapbase_path; print("일반 지도 장면:", mapbase_path())'
  $PY -I osm_extract.py work/osm_features.pkl data/osm/south-korea-latest.osm.pbf data/osm/north-korea-latest.osm.pbf
  $PY -I osm_boundary.py work/territory.pkl data/osm/south-korea-latest.osm.pbf data/osm/north-korea-latest.osm.pbf
  $PY -I korea_mask.py work/korea_mask.npz
  for t in $TILES; do echo ${t%,*} ${t#*,}; done | xargs -P 5 -n 2 sh -c \
    './.venv/bin/python -I terrain.py $0 $1 work/korea_mask.npz work/terrain/$0_$1.json'
  $PY -I land.py
  $PY -I bundle.py
fi
if [[ $HAVE_PREVIEW == 1 ]]; then $PY -I preview.py; fi
