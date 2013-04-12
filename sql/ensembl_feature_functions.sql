select gene.stable_id, description
  from gene
 where description is not null
 union
select transcript.stable_id, description
  from transcript
 where description is not null
 