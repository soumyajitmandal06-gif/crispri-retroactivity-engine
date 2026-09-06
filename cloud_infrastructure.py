import psycopg2
import sys

# ---> REPLACE THIS STRING with the exact URL you copied from Neon <---
NEON_URL = "postgresql://neondb_owner:npg_gDUX5atMH2bz@ep-cool-butterfly-adfbjv0g-pooler.c-2.us-east-1.aws.neon.tech/neondb?sslmode=require&channel_binding=require"

E_COLI_CONSTANTS = [
    ('RNA Polymerase Elongation Rate', 39, 90, 'nt/s', 'BNID_105449', 'E. coli'),
    ('Ribosome Translation Rate', 12, 21, 'aa/s', 'BNID_100059', 'E. coli'),
    ('Total Free Ribosomes', 10000, 60000, 'count', 'BNID_101440', 'E. coli'),
    ('mRNA Half-Life', 2, 5, 'min', 'BNID_100650', 'E. coli')
]

try:
    print("--> Connecting to Neon Cloud Database...")
    conn = psycopg2.connect(NEON_URL)
    cur = conn.cursor()

    print("--> Executing Relational Schema...")
    cur.execute("""
        DROP TABLE IF EXISTS circuit_assemblies, kinetic_constants, genetic_parts CASCADE;

        CREATE TABLE IF NOT EXISTS genetic_parts (
            part_id VARCHAR(100) PRIMARY KEY,
            part_type VARCHAR(50) NOT NULL,
            sequence TEXT,
            source_database VARCHAR(50),
            description TEXT
        );

        CREATE TABLE IF NOT EXISTS kinetic_constants (
            parameter_id SERIAL PRIMARY KEY,
            part_id VARCHAR(100) REFERENCES genetic_parts(part_id),
            parameter_type VARCHAR(50) NOT NULL,
            min_value NUMERIC NOT NULL,
            max_value NUMERIC NOT NULL,
            unit VARCHAR(50) NOT NULL,
            bionumber_id VARCHAR(50) UNIQUE,
            organism_context VARCHAR(100)
        );

        CREATE TABLE IF NOT EXISTS circuit_assemblies (
            assembly_id SERIAL PRIMARY KEY,
            circuit_name VARCHAR(150) NOT NULL,
            promoter_id VARCHAR(100) REFERENCES genetic_parts(part_id),
            rbs_id VARCHAR(100) REFERENCES genetic_parts(part_id),
            cds_id VARCHAR(100) REFERENCES genetic_parts(part_id),
            terminator_id VARCHAR(100) REFERENCES genetic_parts(part_id)
        );
    """)

    print("--> Seeding Baseline Constraints...")
    cur.execute("""
        INSERT INTO genetic_parts (part_id, part_type, source_database, description)
        VALUES ('genome_ecoli_mg1655', 'Host_Genome', 'BioNumbers', 'Baseline E. coli host chassis')
        ON CONFLICT (part_id) DO NOTHING;
    """)

    for param in E_COLI_CONSTANTS:
        cur.execute("""
            INSERT INTO kinetic_constants (part_id, parameter_type, min_value, max_value, unit, bionumber_id, organism_context)
            VALUES ('genome_ecoli_mg1655', %s, %s, %s, %s, %s, %s)
            ON CONFLICT (bionumber_id) DO NOTHING;
        """, param)

    conn.commit()
    print("--> [SUCCESS] Cloud architecture initialized and verified.")

except Exception as e:
    print(f"--> [FATAL ERROR] Pipeline failed: {e}")
    sys.exit(1)

finally:
    if 'conn' in locals() and conn:
        cur.close()
        conn.close()