select gene_stable_id.stable_id, description
  from gene left join gene_stable_id using (gene_id)
 where description is not null
 union
select transcript_stable_id.stable_id, description
  from transcript left join transcript_stable_id using (transcript_id)
 where description is not null
 