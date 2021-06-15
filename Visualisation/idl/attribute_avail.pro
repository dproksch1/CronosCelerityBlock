FUNCTION attribute_avail, group_id, search_name
  is_avail = 0
  n_attr = H5A_Get_Num_Attrs (group_id)
  For iatt = 0, n_attr-1 Do Begin
     attrset_id = H5A_Open_Idx (group_id, iatt)
     attr_name  = H5A_Get_Name (attrset_id)
     If(attr_name eq search_name) Then is_avail = 1
  EndFor
  return,is_avail
END
