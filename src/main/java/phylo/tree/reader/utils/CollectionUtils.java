/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phylo.tree.reader.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 *
 * @author tanxu
 */
public class CollectionUtils {
    
    /**
     * copy all the element in array to an ArrayList object that is not backed by the array;
     * this method is different from Arrays.asList(), which will create List backed by the input array and can not be added more elements;
     * @param <T>
     * @param array
     * @return 
     */
    public static <T> List<T> copyArrayToList(T[] array){
        List<T> ret = new ArrayList<>();
        for(T t:array){
            ret.add(t);
        }
        return ret;
    }
    
    
    
    
    /**
     * check if the given two sets are disjoint with each other (no non-null element in both sets)
     * @param <T>
     * @param set1
     * @param set2
     * @return
     */
    public static <T> boolean setsAreDisjoint(Set<T> set1, Set<T> set2) {
		for(T t: set1) {
			if(t==null) //skip null element
				continue;
			
			if(set2.contains(t))
				return false;
		}
		return true;
	}
	
    /**
     * 
     * @param <T>
     * @param set1
     * @param set2
     * @return
     */
    public static <T> boolean set1IsSuperSetOfSet2(Set<T> set1, Set<T> set2) {
		for(T t:set2) {
			if(!set1.contains(t)) {
				return false;
			}
		}
		return true;
	}
    
}
