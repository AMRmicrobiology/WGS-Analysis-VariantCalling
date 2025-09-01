#!/usr/bin/env bash
set -Eeuo pipefail

# ==========================
# Config
# ==========================
DB_DIR="${DB_DIR:-/kraken2-db}"           # destino por defecto

if [[ "$PWD" != "/" ]]; then
  DB_DIR="${DB_DIR:-$PWD/kraken_db}"
fi

DB_SELECT="${DB_SELECT:-db_16GB}"         # db_full_60GB | db_16GB | /ruta/local | /ruta/db.tar.gz
DB_URL="${DB_URL:-}"                      # opcional: URL directa para override manual
DB_URL_CHECKSUM="${DB_URL_CHECKSUM:-}"    # opcional: URL a .sha256 o .md5

# URLs FIJAS (ajusta cuando haya nuevas releases)
DB_URL_FULL_DEFAULT="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250714.tar.gz"
DB_URL_16_DEFAULT="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20250714.tar.gz"

export KRAKEN2_DEFAULT_DB="$DB_DIR"

# Descarga reanudable en el volumen (no en /tmp) + lock anti-carreras
CACHE_DIR="${CACHE_DIR:-$DB_DIR/.cache}"
TGZ_PATH="$CACHE_DIR/kraken_db.tar.gz"
LOCK_DIR="$DB_DIR/.download.lock"

# ==========================
# Helpers
# ==========================
have_db () {
  [[ -f "$DB_DIR/hash.k2d" && -f "$DB_DIR/opts.k2d" && -f "$DB_DIR/taxo.k2d" ]]
}

download_with_resume () {
  local url="$1"
  echo "Descargando (reanudable) en $TGZ_PATH"
  aria2c -x 8 -s 8 -c -o "$(basename "$TGZ_PATH")" -d "$CACHE_DIR" "$url"
}

verify_checksum_from_url () {
  local tgz="$1" sum_url="$2"
  echo "Verificando checksum: $sum_url"
  if ! curl -fsSL "$sum_url" -o "$CACHE_DIR/sum.chk"; then
    echo "Aviso: no pude obtener checksum, continúo sin verificación."
    return 0
  fi
  case "$sum_url" in
    *.sha256)
      [[ "$(sha256sum "$tgz" | awk '{print $1}')" == "$(awk '{print $1}' "$CACHE_DIR/sum.chk")" ]] \
        || { echo "ERROR: SHA256 no coincide"; exit 6; }
      ;;
    *.md5)
      [[ "$(md5sum "$tgz" | awk '{print $1}')" == "$(awk '{print $1}' "$CACHE_DIR/sum.chk")" ]] \
        || { echo "ERROR: MD5 no coincide"; exit 5; }
      ;;
    *) echo "Aviso: extensión de checksum desconocida, omito verificación." ;;
  esac
}

detect_strip_components () {
  local tgz="$1"
  local first
  first="$(tar -tzf "$tgz" | head -n1 || true)"
  if [[ "$first" == */ ]]; then
    echo 1
  else
    echo 0
  fi
}

verify_db () {
  for f in hash.k2d opts.k2d taxo.k2d; do
    if [[ ! -f "$DB_DIR/$f" ]]; then
      echo "ERROR: falta $f en $DB_DIR"
      echo "Contenido de $DB_DIR:"
      ls -la "$DB_DIR"
      exit 9
    fi
  done
}

extract_db () {
  local tgz="$1"
  echo "Extrayendo en $DB_DIR ..."
  mkdir -p "$DB_DIR"
  local strip
  strip="$(detect_strip_components "$tgz")"
  if [[ "$strip" -eq 1 ]]; then
    tar -C "$DB_DIR" --strip-components=1 -xzf "$tgz"
  else
    tar -C "$DB_DIR" -xzf "$tgz"
  fi
  verify_db
  echo "DB lista en: $DB_DIR"
}

# Selección de fuente según DB_SELECT / DB_URL
prepare_db_from_choice () {
  # 0) Lock anti-carreras
  if ! mkdir "$LOCK_DIR" 2>/dev/null; then
    echo "Otra descarga en curso; esperando a que libere el lock..."
    while [[ -d "$LOCK_DIR" ]]; do sleep 5; done
    return 0  # Al salir del wait, la DB ya debería estar lista
  fi
  trap 'rmdir "$LOCK_DIR" || true' EXIT

  # 1) Override explícito por URL
  if [[ -n "$DB_URL" ]]; then
    download_with_resume "$DB_URL"
    [[ -n "$DB_URL_CHECKSUM" ]] && verify_checksum_from_url "$TGZ_PATH" "$DB_URL_CHECKSUM" || true
    extract_db "$TGZ_PATH"
    return 0
  fi

  # 2) Ruta local a directorio o tarball
  if [[ "$DB_SELECT" == */* ]]; then
    if [[ -d "$DB_SELECT" ]]; then
      echo "Usando DB local en: $DB_SELECT"
      export DB_DIR="$DB_SELECT"
      export KRAKEN2_DEFAULT_DB="$DB_DIR"
      have_db || { echo "ERROR: en '$DB_DIR' no se encuentran hash.k2d/opts.k2d/taxo.k2d"; exit 3; }
      return 0
    elif [[ -f "$DB_SELECT" && "$DB_SELECT" == *.tar.gz ]]; then
      echo "Usando tarball local: $DB_SELECT"
      cp -f "$DB_SELECT" "$TGZ_PATH"
      extract_db "$TGZ_PATH"
      return 0
    else
      echo "ERROR: DB_SELECT parece ruta, pero no existe: $DB_SELECT" >&2
      exit 4
    fi
  fi

  # 3) Presets internos
  local url=""
  case "$DB_SELECT" in
    db_full_60GB) url="$DB_URL_FULL_DEFAULT" ;;
    db_16GB)      url="$DB_URL_16_DEFAULT"   ;;
    *) echo "ERROR: DB_SELECT debe ser 'db_full_60GB', 'db_16GB' o una ruta local"; exit 2 ;;
  esac

  download_with_resume "$url"
  [[ -n "$DB_URL_CHECKSUM" ]] && verify_checksum_from_url "$TGZ_PATH" "$DB_URL_CHECKSUM" || true
  extract_db "$TGZ_PATH"
}

# ==========================
# Entradas de usuario / flujo
# ==========================

# Permite: kraken2-entrypoint.sh prepare-db [preset]
if [[ "${1:-}" == "prepare-db" ]]; then
  if [[ -n "${2:-}" ]]; then DB_SELECT="$2"; fi
  if have_db; then
    echo "DB ya presente en $DB_DIR. Verificando..."
    verify_db
    exit 0
  fi
  echo "No se encontró DB en $DB_DIR. Preparando según DB_SELECT='$DB_SELECT'..."
  prepare_db_from_choice
  exit 0
fi

# Wrapper normal: si no hay DB, prepárala antes de ejecutar kraken2
if ! have_db; then
  echo "No se encontró DB en $DB_DIR. Preparando según DB_SELECT='$DB_SELECT'..."
  prepare_db_from_choice
else
  verify_db
fi

echo "Usando DB en $DB_DIR"
exec kraken2 "$@"
