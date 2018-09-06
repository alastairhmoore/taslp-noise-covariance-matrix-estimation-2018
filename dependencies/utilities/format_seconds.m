function[out_str] = format_seconds(remaining_secs)
minute = 60;
hour = 60*minute;
day = 24*hour;
units = [day,hour,minute];
unit_labels = {'day','hour','min'};
out_str = [];


for iu = 1:length(units)
    [out_str,remaining_secs] = consume_unit(out_str,remaining_secs,units(iu),unit_labels{iu});
end

out_str = [out_str,sprintf(' %2.2f s',remaining_secs)];


function[out_str,remaining_secs] = consume_unit(out_str,remaining_secs,unit,label)
if remaining_secs > unit
    whole_units = floor(remaining_secs/unit);
    if whole_units > 1
        out_str = [out_str, sprintf(' %d %ss',whole_units,label)];
    else
        out_str = [out_str, sprintf(' %d %s',whole_units,label)];
    end
    remaining_secs = remaining_secs - whole_units*unit;
end
