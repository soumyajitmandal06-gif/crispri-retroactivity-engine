import psycopg2
import sys

# ---> PASTE YOUR EXACT NEON URL HERE <---
NEON_URL = "postgresql://neondb_owner:npg_gDUX5atMH2bz@ep-cool-butterfly-adfbjv0g-pooler.c-2.us-east-1.aws.neon.tech/neondb?sslmode=require&channel_binding=require"

def inject_crispri_data():
    try:
        conn = psycopg2.connect(NEON_URL)
        cur = conn.cursor()
        
        print("--> Connecting to Neon Cloud...")
        
        # 1. Inject the new genetic parts (Ignoring if they already exist)
        cur.execute("""
            INSERT INTO genetic_parts (part_id, part_type, description) 
            VALUES 
            ('dcas9_protein', 'Repressor', 'Catalytically dead Cas9'),
            ('sgrna_transcript', 'Guide RNA', 'Single guide RNA targeting primary promoter')
            ON CONFLICT DO NOTHING;
        """)
        
        # 2. Define the exact biophysical constraints for CRISPRi
        crispri_kinetics = [
            ('dcas9_protein', 'Translation Rate', 15.0, 20.0, 'aa/sec'),
            ('dcas9_protein', 'Protein Half-Life', 600, 1200, 'minutes'),
            ('sgrna_transcript', 'Transcription Rate', 0.1, 0.5, 'mRNA/sec'),
            ('sgrna_transcript', 'RNA Half-Life', 2.0, 5.0, 'minutes'),
            ('dcas9_protein', 'Complex Association (kon)', 0.001, 0.01, '1/(nM*s)'),
            ('dcas9_protein', 'Complex Dissociation (koff)', 0.0001, 0.001, '1/s')
        ]
        
        # 3. Inject the kinetics into the database
     # 3. Inject the kinetics into the database
        for part, p_type, min_v, max_v, unit in crispri_kinetics:
            cur.execute("""
                INSERT INTO kinetic_constants (part_id, parameter_type, min_value, max_value, unit) 
                VALUES (%s, %s, %s, %s, %s)
            """, (part, p_type, min_v, max_v, unit))
        
        conn.commit()
        print("--> [SUCCESS] CRISPRi biophysical parameters injected into database.")
        
        cur.close()
        conn.close()
        
    except Exception as e:
        print(f"--> [DATABASE ERROR] {e}")
        sys.exit(1)

if __name__ == "__main__":
    inject_crispri_data()